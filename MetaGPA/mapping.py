from .template import Processor, Settings


class Mapping(Processor):

    fna: str
    control_fq1: str
    control_fq2: str
    case_fq1: str
    case_fq2: str

    index: str

    control_bam: str
    case_bam: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            fna: str,
            control_fq1: str,
            control_fq2: str,
            case_fq1: str,
            case_fq2: str):

        self.fna = fna
        self.control_fq1 = control_fq1
        self.control_fq2 = control_fq2
        self.case_fq1 = case_fq1
        self.case_fq2 = case_fq2

        self.index_genome()
        self.map_reads()

        return self.control_bam, self.case_bam

    def index_genome(self):
        self.index = f'{self.workdir}/bowtie2-index'
        cmd = f'bowtie2-build -f {self.fna} {self.index}'
        self.call(cmd)

    def map_reads(self):
        self.control_bam = f'{self.workdir}/control.bam'
        self.case_bam = f'{self.workdir}/case.bam'
        for fq1, fq2, bam in [
            (self.control_fq1, self.control_fq2, self.control_bam),
            (self.case_fq1, self.case_fq2, self.case_bam)
        ]:
            cmd = [
                f'bowtie2 --threads {self.threads} -x {self.index} -1 {fq1} -2 {fq2} --no-unal | samtools view -@ {self.threads} -bS - | samtools sort -@ {self.threads} - > {bam}',
                f'samtools index {bam}']
            for each in cmd:
                self.call(each)
