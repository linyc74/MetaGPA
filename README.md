# MetaGPA

**Metagenomics Genome and Phenome Association (MetaGPA) Analysis**

### Requirements

Executables in `PATH`:
- SPAdes
- cd-hit
- hmmsearch
- bowtie2
- bedtools 

Python 3:
- numpy
- matplotlib
- seaborn
- scipy
- statsmodels
- pandas
- argparse

### Resource Files

`Pfam-A.hmm.gz` can be downloaded from [Release 33.0](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0)

### Usage

Clone the package:

    git clone https://github.com/linyc74/MetaGPA.git

Run with default settings:

    python MetaGPA \
    -1 <control.1.fq> \
    -2 <control.2.fq> \
    -3 <case.1.fq> \
    -4 <case.2.fq> \
    -p <Pfam-A.hmm.gz>

Other options:

    -l <minimum contig length> [default 0]
    -c <enrichment score cutoff> [default 0.0]

### Input Fastq Files

- Can be compressed `.gz`
- Needs to be trimmed, e.g. by `Cutadapt`

### Output Files

1. `all_contigs.nd.fa` all contigs assembled from case and control samples, redundant contigs removed
2. `all_contigs_pfam_a.gtf` hmm annotation file, containing 9 tab-delimited fields:


    1   seqname     <contig_id>
    2   source      .
    3   feature     CDS
    4   start       <start_bp>
    5   end         <end_bp>
    6   score       .
    7   strand      +/-
    8   frame       0
    9   attribute   name "<query_name>";accession "<query_accession>";description "<query_description>";E_value "<E_value>"

3. `contig_enrichment.csv` table for contig enrichment scores, containing 8 comma-delimited fields:


    1   contig_id      <contig_id>
    2   start          <start_bp>
    3   end            <end_bp>
    4   control_counts <mapped reads number in control sample>
    5   case_counts    <mapped reads number in case sample>
    6   RPKM(Ctrl)     <Reads Per Kb per Million reads in control>
    7   RPKM(Case)     <Reads Per Kb per Million reads in case>
    8   enrichment     <RPKM(Case)/RPKM(Control)>

4. `contig_enrichment.png` stripplot of contig enrichment scores, group by control or case contigs
5. `modified_contigs.fa` modified contigs, i.e. enrichment score >= cutoff
6. `unmodified_contigs.fa` unmodified contigs, i.e. enrichment score < cutoff
7. `modified_contigs_pfam_a.gtf` annotation file of modified contigs
8. `unmodified_contigs_pfam_a.gtf` annotation file of unmodified contigs
9. `pfam_a_statistics.csv` statistics of pfam a domains, containing 6 comma-delimited fields:
   

    1   pfam                      <pfam name>
    2   counts_modified_contigs   <counts on modified contigs>
    3   counts_unmodified_contigs <counts on unmodified contigs>
    4   p_value                   <two sided fisher exact test p value>
    5   corrected_p_value         <bonferroni correction for multitest>
