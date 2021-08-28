import os
from typing import List, Tuple
from .template import Processor, Settings


class Assembly(Processor):

    fq1: str
    fq2: str
    assembly_name: str
    min_contigs_length: int

    assembly_dir: str
    fa_data: List[Tuple[str, str]]
    output_fa: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            fq1: str,
            fq2: str,
            assembly_name: str,
            min_contigs_length: int) -> str:

        self.fq1 = fq1
        self.fq2 = fq2
        self.assembly_name = assembly_name
        self.min_contigs_length = min_contigs_length

        self.set_assembly_dir()
        self.spades_assembly()
        self.filter_by_contig_length()
        self.write_output_fa()

        return self.output_fa

    def set_assembly_dir(self):
        self.assembly_dir = f'{self.workdir}/spades/{self.assembly_name}'
        os.makedirs(self.assembly_dir, exist_ok=True)

    def spades_assembly(self):
        cmd = f'spades.py --meta -1 {self.fq1} -2 {self.fq2} -o {self.assembly_dir} --threads {self.threads} --memory {self.memory}'
        self.call(cmd)

    def filter_by_contig_length(self):
        with open(f'{self.assembly_dir}/contigs.fasta') as f:
            self.fa_data = []
            header = f.readline().strip()
            seq = ''
            while True:
                line = f.readline().strip()
                if line.startswith('>'):
                    self.fa_data.append((header, seq))
                    header = line
                    seq = ''
                elif line == '':
                    self.fa_data.append((header, seq))
                    break
                else:
                    seq = seq + line

    def write_output_fa(self):
        self.output_fa = f'{self.workdir}/{self.assembly_name}.fa'
        with open(self.output_fa, 'w') as f:
            for header, seq in self.fa_data:
                if len(seq) >= self.min_contigs_length:
                    print(f'{header}\n{seq}', file=f)


class CombineAssembly(Processor):

    control_fa: str
    case_fa: str

    output_fa: str

    def main(
            self,
            control_fa: str,
            case_fa: str) -> str:

        self.control_fa = control_fa
        self.case_fa = case_fa

        self.output_fa = f'{self.workdir}/all_contigs.fa'

        with open(self.output_fa, 'w') as writer:
            for name, fa in [
                ('control', self.control_fa),
                ('case', self.case_fa)
            ]:
                with open(fa) as reader:
                    while True:
                        line = reader.readline().strip()
                        if line.startswith('>'):
                            id = line[1:]
                            length = line.split('_')[3]
                            new_header = f'>assembly={name};id={id};length={length}'
                            writer.write(new_header + '\n')
                        elif line == '':
                            break
                        else:
                            writer.write(line + '\n')

        return self.output_fa
