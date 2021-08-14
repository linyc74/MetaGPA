import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multi
from .template import Processor, Settings


class CorrectedFisherTest(Processor):

    merged_pfam_count_csv: str

    output_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, merged_pfam_count_csv: str):
        """
        Step 1: do fisher exact test for each pfam (line)
                  counts_selected_contigs      counts_unselected_contigs
          pfam X           a                             b
          not pfam X       c                             d

          a+c = total pfam counts in modified_contigs
          b+d = total pfam counts in unmodified_contigs

        Step 2: do bonferroni correction
        """

        self.merged_pfam_count_csv = merged_pfam_count_csv

        df = pd.read_csv(self.merged_pfam_count_csv, header=0)
        # total number of counts
        N_modified_contigs = sum(df['counts_modified_contigs'])
        N_unmodified_contigs = sum(df['counts_unmodified_contigs'])
        # test for each row
        for i in range(len(df)):
            a = df.loc[i, 'counts_modified_contigs']
            b = df.loc[i, 'counts_unmodified_contigs']
            c = N_modified_contigs - a
            d = N_unmodified_contigs - b

            table = [[a, b],
                     [c, d]]
            oddsratio, p_value = stats.fisher_exact(table, alternative='two-sided')
            df.loc[i, 'p_value'] = p_value

            colname = 'Enriched/depleted in modified contigs'
            if a / N_modified_contigs > b / N_unmodified_contigs:
                df.loc[i, colname] = 'enriched'
            else:
                df.loc[i, colname] = 'depleted'
        df = df.sort_values('p_value')
        pvals = np.array(df['p_value'], dtype='float64')

        reject, pvals_corrected, alphacSidak, alphacBonf = multi.multipletests(
            pvals,
            alpha=0.05,
            method='b',
            is_sorted=True,
            returnsorted=False
        )
        df['corrected_p_value'] = pvals_corrected
        df = df.sort_values('corrected_p_value')

        self.output_csv = f'{self.outdir}/pfam_a_statistics.csv'
        df.to_csv(self.output_csv, index=False, header=True)
