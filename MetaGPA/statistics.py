import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multi
from .template import Processor, Settings


class CorrectedFisherTest(Processor):

    merged_pfam_count_csv: str

    df: pd.DataFrame
    n_modified_contigs: int
    n_unmodified_contigs: int

    output_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, merged_pfam_count_csv: str):
        """
        [Step 1] Fisher exact test for each pfam (line)
                  counts_selected_contigs      counts_unselected_contigs
          pfam X           a                             b
          not pfam X       c                             d

          a+c = total pfam counts in modified_contigs
          b+d = total pfam counts in unmodified_contigs

        [Step 2] Bonferroni correction
        """

        self.merged_pfam_count_csv = merged_pfam_count_csv

        self.df = self.read_data()
        self.set_number_of_contigs()
        self.fisher_test_for_each_row()
        self.bonferroni_correction()
        self.write_output_csv()

    def read_data(self):
        df = pd.read_csv(self.merged_pfam_count_csv, header=0)
        return df

    def set_number_of_contigs(self):
        self.n_modified_contigs = sum(self.df['counts_modified_contigs'])
        self.n_unmodified_contigs = sum(self.df['counts_unmodified_contigs'])

    def fisher_test_for_each_row(self):
        for i in range(len(self.df)):
            a = self.df.loc[i, 'counts_modified_contigs']
            b = self.df.loc[i, 'counts_unmodified_contigs']
            c = self.n_modified_contigs - a
            d = self.n_unmodified_contigs - b

            table = [[a, b],
                     [c, d]]
            oddsratio, p_value = stats.fisher_exact(table, alternative='two-sided')
            self.df.loc[i, 'p_value'] = p_value

            colname = 'Enriched/depleted in modified contigs'
            if a / self.n_modified_contigs > b / self.n_unmodified_contigs:
                self.df.loc[i, colname] = 'enriched'
            else:
                self.df.loc[i, colname] = 'depleted'

        self.df = self.df.sort_values('p_value')

    def bonferroni_correction(self):
        pvals = np.array(self.df['p_value'], dtype='float64')
        reject, pvals_corrected, alphacSidak, alphacBonf = multi.multipletests(
            pvals,
            alpha=0.05,
            method='b',
            is_sorted=True,
            returnsorted=False
        )
        self.df['corrected_p_value'] = pvals_corrected
        self.df = self.df.sort_values('corrected_p_value')

    def write_output_csv(self):
        self.output_csv = f'{self.outdir}/pfam_a_statistics.csv'
        self.df.to_csv(self.output_csv, index=False, header=True)
