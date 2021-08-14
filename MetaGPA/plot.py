import pandas as pd
import seaborn as sns
from .template import Processor, Settings


class PlotContigEnrichment(Processor):

    enrichment_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, enrichment_csv: str):
        self.enrichment_csv = enrichment_csv

        data = pd.read_csv(self.enrichment_csv,header=0)
        df = pd.DataFrame(data)
        df['assembly'] = ''
        for i in range(len(df)):
            contig_id = df.loc[i,'contig_id']
            assembly = contig_id.split('assembly=')[1].split(';')[0] + '_contigs'
            df.loc[i,'assembly'] = assembly
        ax = sns.stripplot(
            data=df,
            x=df['assembly'],
            y=df['enrichment'],
            order=['control_contigs', 'case_contigs'],
            jitter=True,
            alpha=0.5
        )
        ax.set_ylabel('Enrichment')
        ax.set_xlabel(None)
        ax.set_yscale('log')
        ax.set_ylim((0.0001,10000))
        fig = ax.get_figure()
        fig.set_size_inches(w=5,h=5)
        fig.savefig(f'{self.outdir}/contig_enrichment.png', format='png', dpi=600)
