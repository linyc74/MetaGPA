import pandas as pd
from typing import List, Dict, Tuple
from .template import Processor, Settings


class SelectContigs(Processor):

    fna: str
    gtf: str
    enrichment_csv: str
    enrichment_cut_off: float

    modified_contig_ids: List[str]
    unmodified_contig_ids: List[str]
    fna_dict: Dict[str, str]  # contig_id, seq
    gtf_dict: Dict[str, List[str]]  # contig_id, List[gtf_lines]

    modified_gtf: str
    unmodified_gtf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            fna: str,
            gtf: str,
            enrichment_csv: str,
            enrichment_cut_off: float) -> Tuple[str, str]:

        self.fna = fna
        self.gtf = gtf
        self.enrichment_csv = enrichment_csv
        self.enrichment_cut_off = enrichment_cut_off

        self.set_contig_ids()
        self.set_fna_dict()
        self.write_fnas()
        self.set_gtf_dict()
        self.write_gtfs()

        return self.modified_gtf, self.unmodified_gtf

    def set_contig_ids(self):
        data = pd.read_csv(self.enrichment_csv, header=0)
        df = pd.DataFrame(data)
        self.modified_contig_ids = list(df['contig_id'][df['enrichment'] >= self.enrichment_cut_off])
        self.unmodified_contig_ids = list(df['contig_id'][df['enrichment'] < self.enrichment_cut_off])

    def set_fna_dict(self):
        self.fna_dict = {}
        with open(self.fna) as f:
            header = f.readline().strip()[1:]
            seq = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    self.fna_dict[header] = seq
                    header = line[1:]
                    seq = ''
                else:
                    seq += line
            self.fna_dict[header] = seq

    def write_fnas(self):
        with open(f'{self.outdir}/modified_contigs.fa', 'w') as f:
            for each in self.modified_contig_ids:
                print(f'>{each}', file=f)
                print(f'{self.fna_dict[each]}', file=f)

        with open(f'{self.outdir}/unmodified_contigs.fa', 'w') as q:
            for each in self.unmodified_contig_ids:
                print(f'>{each}', file=q)
                print(f'{self.fna_dict[each]}', file=q)

    def set_gtf_dict(self):
        self.gtf_dict = {}
        with open(self.gtf, 'r') as f:
            for line in f:
                line = line.strip()
                key = line.split('\t')[0]
                self.gtf_dict.setdefault(key, []).append(line)

    def write_gtfs(self):
        self.modified_gtf = f'{self.outdir}/modified_contigs_pfam_a.gtf'
        self.unmodified_gtf = f'{self.outdir}/unmodified_contigs_pfam_a.gtf'

        with open(self.modified_gtf, 'w') as f:
            for each in self.modified_contig_ids:
                if each in self.gtf_dict:
                    for item in self.gtf_dict[each]:
                        print(item, file=f)
        with open(self.unmodified_gtf, 'w') as f:
            for each in self.unmodified_contig_ids:
                if each in self.gtf_dict:
                    for item in self.gtf_dict[each]:
                        print(item, file=f)
