import pandas as pd
from .template import Processor, Settings


class CountNonredundantPfams(Processor):

    input_gtf: str
    output_csv: str

    df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            input_gtf: str,
            output_csv: str):
        """
        count nonredundant pfams for each contigs, if a pfam occurs multiple time in a contig, it will be counted only once
        """

        self.input_gtf = input_gtf
        self.output_csv = output_csv

        self.read_input_gtf()
        self.count_nonredundant()
        self.write_output_csv()

    def read_input_gtf(self):
        data = {
            'contig_id': [],
            'pfam': []
        }
        with open(self.input_gtf) as f:
            for line in f:
                line = line.strip()
                contig_id = line.split('\t')[0]
                attribute = line.split('\t')[8]
                pfam = attribute.split('name "')[1].split('";')[0]
                data['contig_id'].append(contig_id)
                data['pfam'].append(pfam)
        self.df = pd.DataFrame(data)

    def count_nonredundant(self):
        # groupby 'contig_id' and 'pfam' to combine contig with multiple same pfams
        self.df = self.df.groupby(['contig_id', 'pfam']).size().reset_index()
        # groupby 'pfam' to get nonredundent counts
        self.df = self.df.groupby('pfam').size().reset_index().rename(columns={0: 'nonredundant_counts'})
        # sort by descending
        self.df = self.df.sort_values(by='nonredundant_counts', ascending=False)

    def write_output_csv(self):
        self.df.to_csv(self.output_csv, index=False, header=True)


class MergePfamCounts(Processor):

    modified_pfam_count_csv: str
    unmodified_pfam_count_csv: str

    merged_df: pd.DataFrame
    merged_pfam_count_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            modified_pfam_count_csv: str,
            unmodified_pfam_count_csv: str) -> str:

        self.modified_pfam_count_csv = modified_pfam_count_csv
        self.unmodified_pfam_count_csv = unmodified_pfam_count_csv

        self.merge()
        self.write_csv()

        return self.merged_pfam_count_csv

    def merge(self):
        df1 = pd.read_csv(self.modified_pfam_count_csv).rename(
            columns={'nonredundant_counts': 'counts_modified_contigs'})
        df2 = pd.read_csv('unmodified_contigs_pfam_a_counts.csv').rename(
            columns={'nonredundant_counts': 'counts_unmodified_contigs'})
        df = pd.merge(df1, df2, how='outer', on='pfam')
        df = df.fillna(value=0)
        df.astype({'counts_modified_contigs': int, 'counts_unmodified_contigs': int})
        self.merged_df = df

    def write_csv(self):
        self.merged_pfam_count_csv = f'{self.outdir}/merged_pfam_a_counts.csv'
        self.merged_df.to_csv(self.merged_pfam_count_csv, index=False, header=True)
