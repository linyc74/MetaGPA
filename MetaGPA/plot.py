import pandas as pd
import seaborn as sns
from typing import Any
from .template import Processor, Settings


class PlotContigEnrichment(Processor):

    enrichment_csv: str

    df: pd.DataFrame
    ax: Any

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, enrichment_csv: str):
        self.enrichment_csv = enrichment_csv

        self.read_data()
        self.add_assembly_column()
        self.stripplot()
        self.config_and_save_fig()

    def read_data(self):
        data = pd.read_csv(self.enrichment_csv, header=0)
        self.df = pd.DataFrame(data)

    def add_assembly_column(self):
        self.df['assembly'] = ''
        for i in range(len(self.df)):
            contig_id = self.df.loc[i, 'contig_id']
            assembly = contig_id.split('assembly=')[1].split(';')[0] + '_contigs'
            self.df.loc[i, 'assembly'] = assembly

    def stripplot(self):
        self.ax = sns.stripplot(
            data=self.df,
            x=self.df['assembly'],
            y=self.df['enrichment'],
            order=['control_contigs', 'case_contigs'],
            jitter=True,
            alpha=0.5
        )

    def config_and_save_fig(self):
        self.ax.set_ylabel('Enrichment')
        self.ax.set_yscale('log')
        self.ax.set_ylim((0.0001, 10000))
        fig = self.ax.get_figure()
        fig.set_size_inches(w=5, h=5)
        fig.savefig(f'{self.outdir}/contig_enrichment.png', format='png', dpi=600)
