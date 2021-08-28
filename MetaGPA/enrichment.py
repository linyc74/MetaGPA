import pandas as pd
from typing import Dict
from .tools import count_fq_reads
from .template import Processor, Settings


class CalculateEnrichment(Processor):

    fna: str
    control_fq1: str
    case_fq1: str
    control_bam: str
    case_bam: str

    million_reads: Dict[str, float]
    bed: str
    multicov_output: str
    output_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            fna: str,
            control_fq1: str,
            case_fq1: str,
            control_bam: str,
            case_bam: str) -> str:

        self.fna = fna
        self.control_fq1 = control_fq1
        self.case_fq1 = case_fq1
        self.control_bam = control_bam
        self.case_bam = case_bam

        self.set_million_reads()
        self.write_bed()
        self.bedtools_multicov()
        self.write_output_csv()

        return self.output_csv

    def set_million_reads(self):
        self.million_reads = {}
        for key, fq in [
            ('case', self.case_fq1),
            ('control', self.control_fq1)
        ]:
            self.million_reads[key] = count_fq_reads(fq) / 1000000

    def write_bed(self):
        self.bed = f'{self.workdir}/tmp.bed'
        writer = open(self.bed, 'w')
        with open(self.fna) as reader:
            head = reader.readline().strip()[1:]
            seq = ''
            for line in reader:
                line = line.strip()
                if line.startswith('>'):
                    print(f'{head}\t0\t{len(seq)}', file=writer)
                    head = line[1:]
                    seq = ''
                else:
                    seq += line
            print(f'{head}\t0\t{len(seq)}', file=writer)
        writer.close()

    def bedtools_multicov(self):
        self.multicov_output = f'{self.workdir}/tmp.multicov'
        cmd = f'bedtools multicov -bams {self.control_bam} {self.case_bam} -bed {self.bed} > {self.multicov_output}'
        self.call(cmd)

    def write_output_csv(self):
        self.output_csv = f'{self.outdir}/contig_enrichment.csv'
        data = pd.read_csv(
            self.multicov_output,
            sep='\t',
            header=None,
            names=['contig_id', 'start', 'end', 'control_counts', 'case_counts']
        )
        df = pd.DataFrame(data)
        df['RPKM(Ctrl)'] = df['control_counts'] / (df['end'] / 1000) / self.million_reads['control']
        df['RPKM(Case)'] = df['case_counts'] / (df['end'] / 1000) / self.million_reads['case']
        df['enrichment'] = df['RPKM(Case)'] / df['RPKM(Ctrl)']
        df.to_csv(self.output_csv, header=True, index=False)
