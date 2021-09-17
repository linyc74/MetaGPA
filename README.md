# MetaGPA

**MetaGPA (Metagenomics Genome-Phenome Association)** pipeline is used to link genetic information (e.g. genes encoding for enzymes) in metagenome with a dedicated functional phenotype (e.g. DNA cytosine modification).
The main MetaGPA pipeline identifies protein family domains (Pfam, from Pfam-A database, http://pfam.xfam.org/) which are significantly associated with the studied phenotype.
We also provide additional scripts for further analyses of any chosen interested Pfam identified from the main pipeline.
For more details, please read our preprint on https://www.biorxiv.org/content/10.1101/2021.03.23.436658v1.

## Usage

Clone the package:

```bash
git clone https://github.com/linyc74/MetaGPA.git
```

Run with default settings:

```bash
python MetaGPA \
-1 <control.1.fq> \
-2 <control.2.fq> \
-3 <case.1.fq> \
-4 <case.2.fq> \
-p <Pfam-A.hmm.gz>
```

Other options:

```bash
-l <minimum contig length> [default 1000]
-c <enrichment score cutoff> [default 3.0]
```

For more details

```bash
python MetaGPA -h
```

## Requirements

Executables in the environment variable `PATH`:
- SPAdes (3.15.3)
- cd-hit (4.8.1)
- hmmer (3.3.2)
- bowtie2 (2.4.1)
- samtools (1.11)
- bedtools (2.30.0)

Python 3 libraries:
- numpy (1.19.2)
- pandas (1.2.4)
- matplotlib (3.3.2)
- seaborn (0.11.0)
- scipy (1.5.2)
- statsmodels (0.12.0)
- argparse (1.1)

*Versions shown in parentheses have been tested, but not required.*

## Resource Files

`Pfam-A.hmm.gz` can be downloaded from [Release 33.0](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0)

## Input Files

- Can be compressed `.gz`
- Fastq files need to be trimmed, e.g. by `Cutadapt`

## Output Files

1. `all_contigs.nd.fa` all contigs assembled from case and control samples, redundant contigs removed
2. `all_contigs_pfam_a.gtf` hmm annotation file, containing 9 tab-delimited fields:

```
1   seqname     <contig_id>
2   source      .
3   feature     CDS
4   start       <start_bp>
5   end         <end_bp>
6   score       .
7   strand      +/-
8   frame       0
9   attribute   name "<query_name>";accession "<query_accession>";description "<query_description>";E_value "<E_value>"
```

3. `contig_enrichment.csv` table for contig enrichment scores, containing 8 comma-delimited fields:

```
1   contig_id      <contig_id>
2   start          <start_bp>
3   end            <end_bp>
4   control_counts <mapped reads number in control sample>
5   case_counts    <mapped reads number in case sample>
6   RPKM(Ctrl)     <Reads Per Kb per Million reads in control>
7   RPKM(Case)     <Reads Per Kb per Million reads in case>
8   enrichment     <RPKM(Case)/RPKM(Control)>
```

4. `contig_enrichment.png` stripplot of contig enrichment scores, group by control or case contigs
5. `modified_contigs.fa` modified contigs, i.e. enrichment score >= cutoff
6. `unmodified_contigs.fa` unmodified contigs, i.e. enrichment score < cutoff
7. `modified_contigs_pfam_a.gtf` annotation file of modified contigs
8. `unmodified_contigs_pfam_a.gtf` annotation file of unmodified contigs
9. `pfam_a_statistics.csv` statistics of pfam a domains, containing 6 comma-delimited fields:
   
```
1   pfam                      <pfam name>
2   counts_modified_contigs   <counts on modified contigs>
3   counts_unmodified_contigs <counts on unmodified contigs>
4   p_value                   <two sided fisher exact test p value>
5   Enriched/depleted         <enriched/depleted in modified contigs>
6   corrected_p_value         <bonferroni correction for multitest>
```

## Docker

Pull Docker image

```bash
docker pull linyc74/metagpa:latest
```

Run container

```bash
docker run \
--volume "/PATH/TO/OUTDIR":"/PATH/TO/OUTDIR" \
--volume "/PATH/TO/INDIR":"/PATH/TO/INDIR" \
linyc74/metagpa:latest \
--control-fq-1 /PATH/TO/INDIR/control.1.fq.gz \
--control-fq-2 /PATH/TO/INDIR/control.2.fq.gz \
--case-fq-1 /PATH/TO/INDIR/case.1.fq.gz \
--case-fq-2 /PATH/TO/INDIR/case.2.fq.gz \
--pfam-hmm /PATH/TO/INDIR/Pfam-A.hmm.gz \
--min-contig-length 1000 \
--enrichment-cutoff 3.0 \
--threads 32 \
--memory 128 \
--outdir /PATH/TO/OUTDIR
```

## Test

To test the MetaGPA main pipeline, a small dataset is provided in the `test` folder, including the following files:

- test_control.1.fq.gz
- test_control.2.fq.gz
- test_case.1.fq.gz
- test_case.2.fq.gz
- Pfam-A-sampled.hmm.gz
