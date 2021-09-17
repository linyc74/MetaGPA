#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from subprocess import check_call

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
def grep_pfam_gtf(input_gtf):
    cmd = f'grep -i "{args.pfam}" {input_gtf} > {input_gtf}.tmp'
    check_call(cmd, shell=True)

# get dna sequence and translate to protein .faa 
def get_faa(gtf,fasta):
    writer = open(f'{gtf}.tmp.faa','w')    
    # read fasta into a dict{'assembly=control;id=NODE_440_length_49190_cov_245.080228;length=49190':'ATCG....',...}
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
    with open(f'{gtf}.tmp','r') as f:
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
            new_header = f'>{contig_id};pfam={args.pfam};start={start};end={end};strand={strand}'
            print(new_header, file=writer)
            print(aa_seq, file=writer)
    writer.close()

# rename header
def rename(faa):
    if 'unmodified' in faa:
        prefix = 'UM'
    else:
        prefix = 'M'
    index_file_name = faa + '.index'
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
def muscle():
    cmd = f'muscle -in {args.pfam}.faa -out {args.pfam}_aligned.afa'
    check_call(cmd, shell=True)

# Raxml
def RAXML():
    cmd = f'raxmlHPC -f a -# autoMRE -p 1237 -x 1237 -m PROTGAMMAAUTO -s {args.pfam}_aligned.afa -n {args.pfam}.tree'
    check_call(cmd, shell=True)

def clear():
    cmd = f'rm {args.pfam}_aligned* *tmp* {args.pfam}.faa'
    check_call(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pfam', help='pfam name', dest='pfam')
    args = parser.parse_args()

    with open('all_contigs_pfam_a.gtf','r') as f:
        text = f.read()

    if not args.pfam in text:
        print('Invalid pfam name')
        quit()

files = [('modified_contigs_pfam_a.gtf','modified_contigs.fa'), 
         ('unmodified_contigs_pfam_a.gtf','unmodified_contigs.fa')]

for gtf, fasta in files:
    grep_pfam_gtf(gtf)
    get_faa(gtf,fasta)
    rename(faa=f'{gtf}.tmp.faa')
    
cmd = f'cat rename_modified*.faa rename_unmodified*.faa > {args.pfam}.faa'
check_call(cmd, shell=True)

cmd = f'cat *.index > {args.pfam}.index'
check_call(cmd, shell=True)
    
muscle()
RAXML()
clear()

# iTOL visualizes RAxML_bipartitions file