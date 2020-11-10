from subprocess import check_call
from os import chdir

def rev_comp(dna):
    D = {'A':'T','T':'A','C':'G','G':'C'}
    rc_dna = ''
    for i in range(len(dna)):
        nt = D.get(dna[len(dna)-1])
        rc_dna += nt
        dna = dna[:len(dna)-1]
    return rc_dna

def translate(dna):
    """
    Translate a DNA sequence.
    If the DNA length is not a multiple of 3, then leave the last 1 or 2 bases untranslated.

    Args:
        dna: str,
            the DNA sequence

    Returns: str,
        the translated amino acid sequence
    """
    codon = {
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'CTT': 'L',
    'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L',
    'TTG': 'L', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
    'GTG': 'V', 'TTT': 'F', 'TTC': 'F', 'ATG': 'M',
    'TGT': 'C', 'TGC': 'C', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GGT': 'G', 'GGC': 'G',
    'GGA': 'G', 'GGG': 'G', 'CCT': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TGG': 'W', 'CAA': 'Q',
    'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'CAT': 'H',
    'CAC': 'H', 'GAA': 'E', 'GAG': 'E', 'GAT': 'D',
    'GAC': 'D', 'AAA': 'K', 'AAG': 'K', 'CGT': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R',
    'AGG': 'R', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

    dna = dna.upper()
    peptide = []
    for i in range(int(len(dna)/3)):
        aa = codon.get(dna[(i*3):(i*3+3)], 'X')
        peptide.append(aa)
    return ''.join(peptide)

# grep pfam of interest
def grep_pfam_gtf(input_gtf, subset_gtf, pfam):
    cmd = f'grep -i "{pfam}" {input_gtf} > {subset_gtf}'
    check_call(cmd, shell=True)

# get dna sequence from .fa and translate to protein .faa (frame was already considered when gtf was writen(start and end of nucleotide))
def get_faa(subset_gtf, fasta, pfam, output_faa):
    # output e.g. Carbam_trans_C_selected_contigs.faa
    writer = open(output_faa,'w')    
    # read fasta into a dict{'source=19_0131_sewage;assembly=control;id=NODE_440_length_49190_cov_245.080228;length=49190':'ATCG....',...}
    dict = {}
    with open(fasta,'r') as f:
        header = f.readline().strip()[1:]
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                dict[header]= seq
                header = line[1:]
                seq = ''
            else:
                seq += line
        dict[header] = seq
    # find the nucleotide sequence of each pfam record and translate
    with open(subset_gtf,'r') as f:
        for line in f:
            line=line.strip()
            contig_id = line.split('\t')[0]
            start = int(line.split('\t')[3])
            end = int(line.split('\t')[4])
            strand = line.split('\t')[6]
            ref_seq = dict[contig_id]       # get the ref sequence
            # get nucleotide sequence, note strand conditions
            if strand == '+':                     
                nt_seq = ref_seq[start-1:end]
            else:
                nt_seq = rev_comp(ref_seq[start-1:end])
            aa_seq = translate(nt_seq)
            # write .faa record
            new_header = f'>{contig_id};pfam={pfam};start={start};end={end};strand={strand}'
            print(new_header, file=writer)
            print(aa_seq, file=writer)
    writer.close()

# rename header
def rename(faa):
    if 'unselected' in faa:
        prefix = 'US'
    else:
        prefix = 'S'
    index_file_name = faa.split('.faa')[0] + '.index'
    index = open(index_file_name,'w')
    writer = open(f'rename_{faa}','w')
    with open(faa,'r') as f:
        i = 1
        header = f.readline().strip()
        seq=''
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                new_header = f'>{prefix}_{i}'
                print(f'{new_header[1:]}\t{header[1:]}',file=index)
                i += 1
                print(new_header, file=writer)
                print(seq, file = writer)
                header = line
                seq = ''
            else:
                seq += line
        new_header = f'>{prefix}_{i}'
        print(f'{new_header[1:]}\t{header[1:]}',file=index)
        print(new_header, file=writer)
        print(seq, file = writer)
    index.close()
    writer.close()
    
# MUSCLE alignment
def muscle(input_faa, aligned_fasta):
    cmd = f'muscle -in {input_faa} -out {aligned_fasta}'
    check_call(cmd, shell=True)

# Raxml
def RAXML(input, output):
    cmd = f'raxmlHPC -f a -# 100 -p 1237 -x 1237 -m PROTGAMMAAUTO -s {input} -n {output}'
    check_call(cmd, shell=True)

chdir('/Users/wyang/Desktop/phage_CT_project/dry_lab/042320_remove_redundant_contigs/phylogenetic_tree')
pfam = 'Carbam_trans_C'
files = [('selected_contigs_pfam_a_95.gtf', 'selected_contigs_95.fa'),
         ('unselected_contigs_pfam_a_95.gtf', 'unselected_contigs_95.fa')]
for gtf, fasta in files:
    subset_gtf = f'{pfam}_{gtf}'
    grep_pfam_gtf(gtf, subset_gtf, pfam)
    protein_faa = subset_gtf.split('_pfam_a_95.gtf')[0] +'.faa'
    get_faa(subset_gtf, fasta, pfam, protein_faa)
    rename(protein_faa)
    
cmd = f'cat rename_*_selected_contigs.faa rename_*_unselected_contigs.faa > {pfam}.faa'
check_call(cmd, shell=True)

cmd = f'cat *.index > {pfam}.index'
check_call(cmd, shell=True)
    
muscle(f'{pfam}.faa', f'{pfam}_aligned.afa')
RAXML(f'{pfam}_aligned.afa',f'{pfam}.tree')

# iTOL visualize RAxML_bipartitions file



