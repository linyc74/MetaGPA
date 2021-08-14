from os.path import exists
from subprocess import check_output


def translate(dna):
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


def rev_comp(dna):
    D = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rc_dna = ''
    for i in range(len(dna)):
        nt = D.get(dna[len(dna)-1])
        rc_dna += nt
        dna = dna[:len(dna)-1]
    return rc_dna


def count_fq_reads(fq: str):
    if '.gz' in fq:
        cmd = f'echo $(zcat {fq} | wc -l)/4 | bc'
        return int(check_output(cmd, shell=True).split()[0])
    else:
        cmd = f'echo $(cat {fq} | wc -l)/4 | bc'
        return int(check_output(cmd, shell=True).split()[0])


def get_temp_path(prefix: str = 'temp', suffix: str = '') -> str:
    i = 0
    while True:
        path = f'{prefix}_{i:06d}{suffix}'
        if not exists(path):
            break
        i += 1
    return path
