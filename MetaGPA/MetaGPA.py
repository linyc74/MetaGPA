import os
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi
from scipy import stats
from scipy.stats import gaussian_kde
from typing import Dict, List
from .template import Processor, Settings


class MetaGPA(Processor):

    GROUPS = ['filter', 'lambda']
    SAMPLES = ['control', 'dnasei']
    CURRENT_DIR = '/mnt/home/ettwiller/wyang/data/101220_DNaseI'
    KMER_SIZES = [14, 16]
    RANDOM_SAMPLE_FOLD = 1000

    million_reads: Dict[str, float]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self):
        for group in self.GROUPS:
            for each in self.SAMPLES:
                self.prepare_file(
                    read1=f'{each}_{group}_reads.1',
                    read2=f'{each}_{group}_reads.2',
                    prefix=f'{each}_{group}')

                for size in self.KMER_SIZES:
                    self.kmer_count(kmer=size, jf_name=f'{each}_{group}')

            for size in self.KMER_SIZES:
                self.kmer_count_to_tsv(
                    control=f'control_{group}.fa',
                    dnasei=f'dnasei_{group}.fa',
                    kmer=size)

        self.million_reads = {}
        for group in self.GROUPS:
            for sample in self.SAMPLES:
                name = f'{sample}_{group}'
                self.million_reads[name] = self.count_reads_per_million(f'{sample}_{group}.fq')

        for kmer in self.KMER_SIZES:
            os.chdir(f'{self.CURRENT_DIR}/kmer_count/{kmer}_mer')

            self.random_sample(
                input='control_filter_vs_dnasei_filter.tsv',
                reduced=f'filter_reduced{self.RANDOM_SAMPLE_FOLD}.tsv',
                fold=self.RANDOM_SAMPLE_FOLD)

            self.density_plot(
                png_name=f'filter_{kmer}_mer_1000.png',
                reduced=f'filter_reduced{self.RANDOM_SAMPLE_FOLD}.tsv',
                group='filter',
                min=0.01,
                max=1000)

            self.density_plot(
                png_name=f'lambda_{kmer}_mer.png',
                reduced=f'control_lambda_vs_dnasei_lambda.tsv',
                group='lambda',
                min=0.01,
                max=100000)

        for sample in self.SAMPLES:
            self.assembly(
                fq1=f'{sample}_filter_reads.1.fq',
                fq2=f'{sample}_filter_reads.2.fq',
                assembly_name=f'{sample}_filter')
            self.edit_assembly(file=f'{sample}_filter_contigs')

        self.combine_assembly()
        self.remove_redundant_contigs()
        self.annotation()
        self.map_reads()
        self.calculate_enrichment()
        self.plot()
        self.select_contigs(threshold=3)

        for sample in self.SAMPLES:
            self.count_nonredundant_pfams(
                input=f'{sample}_contigs_pfam_a.gtf',
                output=f'{sample}_contigs_pfam_a_counts.csv')

        self.merge(
            input_select='select_contigs_pfam_a_counts.csv',
            input_unselect='unselect_contigs_pfam_a_counts.csv',
            output='merged_pfam_a_counts.csv')

        self.corrected_fisher_test(
            input='merged_pfam_a_counts.csv',
            output=f'{self.CURRENT_DIR}/pfam_a_enrichment/pfam_a_enrichment_corrected.csv')

        self.clean()

    def prepare_file(self, read1: str, read2: str, prefix: str):
        os.chdir(f'{self.CURRENT_DIR}/filter')
        cmd = f'cat {read1} {read2} > {prefix}.fq'
        self.call(cmd)

    def kmer_count(self, kmer: int, jf_name: str):
        cmd = f'jellyfish count {prefix}.fq -m {kmer} -s 1G -t 32 -C -o {self.CURRENT_DIR}/kmer_count/{kmer}_mer/{jf_name}.jf'
        self.call(cmd)

        cmd = f'jellyfish dump -L 10 {self.CURRENT_DIR}/kmer_count/{kmer}_mer/{jf_name}.jf > {self.CURRENT_DIR}/kmer_count/{kmer}_mer/{jf_name}.fa'
        self.call(cmd)

    def kmer_count_to_tsv(self, control: str, dnasei: str, kmer: int):
        os.chdir(f'{self.CURRENT_DIR}/kmer_count/{kmer}_mer')
        cmd = f'python {self.CURRENT_DIR}/kmer_count_to_tsv.py {control} {dnasei}'
        self.call(cmd)

    def count_reads_per_million(self, file: str) -> float:
        os.chdir(f'{self.CURRENT_DIR}/filter')
        line_count = 0
        with open(file, 'r') as f:
            for line in f:
                line_count += 1
        reads_count = line_count / 4
        return reads_count / 1000000

    def random_sample(self, input: str, reduced: str, fold: int):
        with open(input, 'r') as f:
            head = f.readline().strip()
            lines = f.read().split('\n')[:-1]

        length = len(lines)
        new = random.sample(lines, int(length / fold))

        with open(reduced, 'w') as writer:
            print(head, file=writer)
            for each in new:
                print(each, file=writer)

    def density_plot(
            self,
            png_name: str,
            reduced: str,
            group: str,
            min: float,
            max: float):

        sequence, count1, count2 = self.read_tsv_into_arrays(reduced)

        # kmer counts per million reads
        x = count1 / self.million_reads[f'control_{group}']
        y = count2 / self.million_reads[f'dnasei_{group}']

        # show data points on axes (replace 0 to 0.03)
        x[x == 0] = 0.03
        y[y == 0] = 0.03

        # calculate point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)

        # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

        # plot
        plt.figure(figsize=(5, 5), dpi=600)
        pcm = plt.scatter(x, y, c=z, cmap='coolwarm', s=1, edgecolor='')
        plt.xlabel('Control kmer counts per million reads')
        plt.ylabel('DNaseI kmer counts per million reads')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(min, max)
        plt.ylim(min, max)
        plt.colorbar(pcm)
        plt.savefig(png_name, format='png', dpi=600)

    def read_tsv_into_arrays(self, reduced: str):
        data = []
        with open(reduced, 'r') as f:
            header = f.readline()
            for line in f:
                data.append(tuple(line.strip().split('\t')))  # data is list containing tuples (kmer, count1, count2) as items
        kmer, count1, count2 = zip(*data)      # unzip into three tuples (by column)
        return kmer, np.array(count1, dtype='int32'), np.array(count2, dtype='int32')

    def assembly(self, fq1: str, fq2: str, assembly_name: str):
        os.chdir(f'{self.CURRENT_DIR}/filter')
        cmd = f'spades.py --meta -1 {fq1} -2 {fq2} -o {self.CURRENT_DIR}/assembly_output/{assembly_name}_contigs --threads {self.threads} --memory 350'
        self.call(cmd)

    def edit_assembly(self, file: str):
        os.chdir(f'{self.CURRENT_DIR}/assembly_output/{file}')
        # read original .fasta into a list of tuples [(header, seq), ...]
        with open('contigs.fasta', 'r') as f:
            temp_list = []
            header = f.readline().strip()
            seq = ''
            while True:
                line = f.readline().strip()
                if line.startswith('>'):
                    temp_list.append((header, seq))
                    header = line
                    seq = ''
                elif line == '':
                    temp_list.append((header, seq))
                    break
                else:
                    seq = seq + line

        min_contigs_length = 1000

        # write new .fasta containing contigs >= min_contigs_length
        os.chdir(f'{self.CURRENT_DIR}/assembly_output/output_contigs')
        file_name = each.split('/')[-1]
        with open(f'{file_name}.fa', 'w') as f:
            for header, seq in temp_list:
                if len(seq) >= min_contigs_length:
                    print(f'{header}\n{seq}', file=f)
                else:
                    continue

    def combine_assembly(self):
        os.chdir(f'{self.CURRENT_DIR}/assembly_output/output_contigs')

        control, dnasei = 'control_filter_contigs.fa', 'dnasei_filter_contigs.fa'

        with open('all_filter_contigs.fa', 'w') as f:
            with open(dnasei, 'r') as f1:
                assembly = 'dnasei'
                while True:
                    line = f1.readline().strip()
                    if line.startswith('>'):
                        id = line[1:]
                        length = line.split('_')[3]
                        new_header = f'>assembly={assembly};id={id};length={length}'
                        print(new_header, file = f)
                    elif line == '':
                        break
                    else:
                        print(line, file=f)
                # read the second file
            with open(control, 'r') as f2:
                assembly = 'control'
                while True:
                    line = f2.readline().strip()
                    if line.startswith('>'):
                        id = line[1:]
                        length = line.split('_')[3]
                        new_header = f'>assembly={assembly};id={id};length={length}'
                        print(new_header, file = f)
                    elif line == '':
                        break
                    else:
                        print(line, file=f)

    def remove_redundant_contigs(self):
        cmd = f'''\
cd-hit-est \
-i {self.CURRENT_DIR}/assembly_output/output_contigs/all_filter_contigs.fa \
-o {self.CURRENT_DIR}/nonredundant_contigs/nonredundant_all_filter_contigs.fa \
-c 0.95 -n 10 -d 0 -M 0 \
-T {self.threads}'''
        self.call(cmd)

    def annotation(self):
        cmd = f'python {self.CURRENT_DIR}/translation.py'
        self.call(cmd)

        cmd = f'''\
hmmsearch \
--cpu {self.threads} \
-o {self.CURRENT_DIR}/pfam_a_annotation/nonredundant_all_filter_contigs.out \
/mnt/home/ettwiller/wyang/data/032420_reanalysis/pfam_a_annotation/Pfam-A.hmm \
{self.CURRENT_DIR}/nonredundant_contigs/nonredundant_all_filter_contigs_translated.fa'''
        self.call(cmd)

        cmd = f'python {self.CURRENT_DIR}/parse_hmm_result_pfam_a.py'
        self.call(cmd)

    def map_reads(self):
        MapReads(self.settings).main(
            current_dir=self.CURRENT_DIR,
            samples=self.SAMPLES)

    def calculate_enrichment(self):
        # generate bed files
        os.chdir(f'{self.CURRENT_DIR}/nonredundant_contigs/')
        bed = open('temp.bed', 'w')
        with open('nonredundant_all_filter_contigs.fa', 'r') as f:
            head = f.readline().strip()[1:]
            seq = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    print(f'{head}\t0\t{len(seq)}', file=bed)
                    head = line[1:]
                    seq = ''
                else:
                    seq += line
            print(f'{head}\t0\t{len(seq)}', file=bed)   # don't forget the last contig
        bed.close()
        # count reads with bedtools multicov
        bam1 = f'{self.CURRENT_DIR}/mapping/control_map_all_filter_contigs.bam'
        bam2 = f'{self.CURRENT_DIR}/mapping/dnasei_map_all_filter_contigs.bam'
        cmd = f'bedtools multicov -bams {bam1} {bam2} -bed temp.bed > multicov'
        self.call(cmd)
        # read into pandas, add header
        data = pd.read_csv('multicov', sep='\t', header=None, names=['contig_id','start','end','control_counts','dnasei_counts'])
        df = pd.DataFrame(data)
        # compute RPKM = Reads Per Kb per Million reads
        df['RPKM(C)'] = df['control_counts'] / (df['end']/1000) / self.million_reads['control']
        df['RPKM(D)'] = df['dnasei_counts'] / (df['end']/1000) / self.million_reads['dnasei']
        df['enrichment'] = df['RPKM(D)'] / df['RPKM(C)']
        # write to csv
        df.to_csv(f'{self.CURRENT_DIR}/select_contigs/contig_enrichment.csv', header=True, index=False)

    def plot(self): # plot contig enrichment
        data = pd.read_csv(f'{self.CURRENT_DIR}/select_contigs/contig_enrichment.csv',header=0)
        df = pd.DataFrame(data)
        df['assembly']=''
        for i in range(len(df)):
            contig_id = df.loc[i,'contig_id']
            assembly = contig_id.split('assembly=')[1].split(';')[0] + '_contigs'
            df.loc[i,'assembly'] = assembly
        ax = sns.stripplot(data=df, x=df['assembly'], y=df['enrichment'], order=['control_contigs','dnasei_contigs'], jitter=True, alpha=1)
        ax.set_ylabel('Enrichment')
        ax.set_xlabel(None)
        ax.set_yscale('log')
        ax.set_ylim((0.0001,10000))
        fig = ax.get_figure()
        fig.set_size_inches(w=5,h=5)
        fig.savefig(f'{self.CURRENT_DIR}/select_contigs/enrichment_log.png', format='png', dpi=600)

    def select_contigs(self, threshold: float):
        data = pd.read_csv(f'{self.CURRENT_DIR}/select_contigs/contig_enrichment.csv',header=0)
        df = pd.DataFrame(data)
        selected_contig_id = list(df['contig_id'][df['enrichment']>=threshold])
        unselected_contig_id = list(df['contig_id'][df['enrichment']<threshold])
        print(f'selected contigs: {len(selected_contig_id)}; unselected contigs: {len(unselected_contig_id)}; total contigs: {len(df)}')
        # generate select_contigs.fa and unselect_contigs.fa
        fasta={}     # read all_filter_contigs.fa and create a dict with key,value = contig_id, seq
        with open(f'{self.CURRENT_DIR}/nonredundant_contigs/nonredundant_all_filter_contigs.fa', 'r') as f:
            header = f.readline().strip()[1:]
            seq = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    fasta[header]=seq
                    header = line[1:]
                    seq = ''
                else:
                    seq += line
            fasta[header]=seq
        # generate selected_contigs.fa by searching each in selected_contig_id
        with open(f'{self.CURRENT_DIR}/select_contigs/select_contigs.fa','w') as f:
            for each in selected_contig_id:
                print(f'>{each}',file = f)
                print(f'{fasta[each]}', file=f)
        # generate unselected_contigs.fa by searching each in unselected_contig_id
        with open(f'{self.CURRENT_DIR}/select_contigs/unselect_contigs.fa','w') as q:
            for each in unselected_contig_id:
                print(f'>{each}',file = q)
                print(f'{fasta[each]}', file=q)
        # generate select_contigs_pfam_a.gtf and unselect_contigs_pfam_a.gtf
        gtf = {}  # read the gtf file into a dictionary with key,value = contig_id, line  # note here each key may have multiple values
        with open(f'{self.CURRENT_DIR}/pfam_a_annotation/nonredundant_all_filter_contigs_pfam_a.gtf','r') as f:
            for line in f:
                line = line.strip()
                key = line.split('\t')[0]
                gtf.setdefault(key, []).append(line)
        # generate selected_contigs_pfam.gtf by search each in selected_contig_id
        with open(f'{self.CURRENT_DIR}/pfam_a_annotation/select_contigs_pfam_a.gtf','w') as f:
            for each in selected_contig_id:
                if each in gtf:
                    for item in gtf[each]:
                        print(item, file=f)
        # generate selected_contigs_pfam.gtf by search each in selected_contig_id
        with open(f'{self.CURRENT_DIR}/pfam_a_annotation/unselect_contigs_pfam_a.gtf','w') as f:
            for each in unselected_contig_id:
                if each in gtf:
                    for item in gtf[each]:
                        print(item, file=f)

    def count_nonredundant_pfams(self, input: str, output: str):
        '''
        count nonredundant pfams for each contigs, if a pfam occurs multiple time in a contig, it will be counted only once
        '''
        # read input data in a dataframe
        os.chdir(f'{self.CURRENT_DIR}/pfam_a_annotation')
        data = {'contig_id': [], 'pfam': []}
        with open(input, 'r') as f:
            for line in f:
                line = line.strip()
                contig_id = line.split('\t')[0]
                attribute = line.split('\t')[8]
                pfam = attribute.split('name "')[1].split('";')[0]
                data['contig_id'].append(contig_id)
                data['pfam'].append(pfam)
        df = pd.DataFrame(data)
        # groupby 'contig_id' and 'pfam' to combine contig with multiple same pfams
        df = df.groupby(['contig_id', 'pfam']).size().reset_index()
        # groupby 'pfam' to get nonredundent counts
        df = df.groupby('pfam').size().reset_index().rename(columns={0: 'nonredundant_counts'})
        # sort by descending and save as csv
        df = df.sort_values(by='nonredundant_counts', ascending=False)
        df.to_csv(output, index=False, header=True)

    def merge(self, input_select: str, input_unselect: str, output: str):
        """
        merge counts in selected_contigs and unselected_contigs

        Args:
        two input files: counts.csv file     e.g.  pfam nonredundant_counts
        output: csv file                     e.g.  pfam counts_selected_contigs counts_unselected_contigs
        """
        os.chdir(f'{self.CURRENT_DIR}/pfam_a_annotation')
        df1 = pd.read_csv(input_select).rename(columns={'nonredundant_counts': 'counts_selected_contigs'})
        df2 = pd.read_csv(input_unselect).rename(columns={'nonredundant_counts': 'counts_unselected_contigs'})
        df = pd.merge(df1, df2, how='outer',
                      on='pfam')  # merge by 'pfam', outer means use union of keys, save all data in two df
        df = df.fillna(value=0)
        df.astype({'counts_selected_contigs': int, 'counts_unselected_contigs': int})
        df.to_csv(output, index=False, header=True)

    def corrected_fisher_test(self, input: str, output: str):
        """
        Fisher exact test for each pfam domain and bonferroni correction for multiple testing

        do fisher exact test for each pfam (line)
                  counts_selected_contigs      counts_unselected_contigs
          pfam X           a                             b
          not pfam X       c                             d

          a+c = total pfam counts in selected_contigs
          b+d = total pfam counts in unselected_contigs

          do bonferroni correction
        output: add p_value and corrected p_value column
        """
        os.chdir(f'{self.CURRENT_DIR}/pfam_a_annotation')
        df = pd.read_csv(input, header=0)
        # total number of counts
        N_selected_contigs = sum(df['counts_selected_contigs'])
        N_unselected_contigs = sum(df['counts_unselected_contigs'])
        # test for each row
        for i in range(len(df)):
            a = df.loc[i,'counts_selected_contigs']
            b = df.loc[i,'counts_unselected_contigs']
            c = N_selected_contigs - a
            d = N_unselected_contigs - b

            table = [[a,b],
                     [c,d]]
            oddsratio, p_value = stats.fisher_exact(table, alternative='two-sided')
            df.loc[i, 'p_value'] = p_value

            colname = 'Enriched/depleted in selected contigs'
            if a/N_selected_contigs > b/N_unselected_contigs:
                df.loc[i,colname] = 'enriched'
            else:
                df.loc[i,colname] = 'depleted'
        df = df.sort_values('p_value')
        pvals = np.array(df['p_value'], dtype='float64')
        reject, pvals_corrected, alphacSidak, alphacBonf = multi.multipletests(pvals, alpha=0.05, method='b', is_sorted=True, returnsorted=False)
        df['corrected_p_value'] = pvals_corrected
        df = df.sort_values('corrected_p_value')
        df.to_csv(output, index=False, header=True)

    def clean(self):
        os.chdir(f'{self.CURRENT_DIR}/pfam_a_annotation')
        self.call('rm *_counts.csv')


class MapReads(Processor):

    current_dir: str
    samples: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            current_dir: str,
            samples: List[str]):

        self.current_dir = current_dir
        self.samples = samples

        self.build_index()
        for sample in self.samples:
            self.mapping(sample)
            self.sam_to_bam(sample)
            self.sort_bam(sample)
            self.index_bam(sample)
        self.clean_up()

    def build_index(self):
        cmd = f'''\
bowtie2-build \
-f {self.current_dir}/nonredundant_contigs/nonredundant_all_filter_contigs.fa \
1019_filter'''
        self.call(cmd)

    def mapping(self, sample: str):
        cmd = f'''\
bowtie2 \
-x 1019_filter \
-1 {self.current_dir}/filter/{sample}_filter_reads.1.fq \
-2 {self.current_dir}/filter/{sample}_filter_reads.2.fq \
-S {self.current_dir}/mapping/{sample}_map_all_filter_contigs.sam \
--no-unal'''
        self.call(cmd)

    def sam_to_bam(self, sample: str):
        cmd = f'''\
samtools view -S -b {self.current_dir}/mapping/{sample}_map_all_filter_contigs.sam \
> {sample}_temp'''
        self.call(cmd)

    def index_bam(self, sample: str):
        cmd = f'''samtools index {self.current_dir}/mapping/{sample}_map_all_filter_contigs.bam'''
        self.call(cmd)

    def sort_bam(self, sample: str):
        cmd = f'''\
samtools sort {sample}_temp \
> {self.current_dir}/mapping/{sample}_map_all_filter_contigs.bam'''
        self.call(cmd)

    def clean_up(self):
        cmd = 'rm control_temp dnasei_temp *.bt2'
        self.call(cmd)
