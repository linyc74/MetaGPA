from .template import Processor, Settings
from .tools import translate, rev_comp


class Annotation(Processor):

    fna: str
    pfama: str

    faa: str
    hmm_txt: str
    output_gtf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str, pfama: str) -> str:

        self.fna = fna
        self.pfama = pfama

        self.translate_dna_to_protein()
        self.hmmsearch()
        self.parse_hmm()

        return self.output_gtf

    def translate_dna_to_protein(self):
        self.faa = f'{self.outdir}/translated.fa'
        with open(self.faa, 'w') as writer:
            with open(self.fna) as reader:
                header = reader.readline().strip()
                seq = ''
                protein_seq = ''
                while True:
                    line = reader.readline()
                    if line.strip().startswith('>'):
                        for frame in (1, 2, 3):
                            new_header = f'{header};frame={frame}'
                            protein_seq = translate(seq[frame - 1:])
                            print(new_header, file=writer)
                            print(protein_seq, file=writer)
                        for frame in (-1, -2, -3):
                            new_header = f'{header};frame={frame}'
                            rc_seq = rev_comp(seq)
                            protein_seq = translate(rc_seq[abs(frame) - 1:])
                            print(new_header, file=writer)
                            print(protein_seq, file=writer)
                        header = line.strip()
                        seq = ''
                    elif line == '':
                        for frame in (1, 2, 3):
                            new_header = f'{header};frame={frame}'
                            protein_seq = translate(seq[frame - 1:])
                            print(new_header, file=writer)
                            print(protein_seq, file=writer)
                        for frame in (-1, -2, -3):
                            new_header = f'{header};frame={frame}'
                            rc_seq = rev_comp(seq)
                            protein_seq = translate(rc_seq[abs(frame) - 1:])
                            print(new_header, file=writer)
                            print(protein_seq, file=writer)
                        break
                    else:
                        seq = seq + line.strip()

    def hmmsearch(self):
        self.hmm_txt = f'{self.workdir}/all_contigs.out'
        cmd = f'hmmsearch --cpu {self.threads} -o {self.hmm_txt} {self.pfama} {self.faa}'
        self.call(cmd)

    def parse_hmm(self):
        ''' parse the result by hmmer and create a gtf file with 9 tab-delimited fields:
                    1   seqname     <contig_id>
                    2   source      .
                    3   feature     CDS
                    4   start       <start_bp>
                    5   end         <end_bp>
                    6   score       .
                    7   strand      +/-
                    8   frame       0
                    9   attribute   name "<query_name>";accession "<query_accession>";description "<query_description>";E_value "<E_value>"
        '''

        with open(self.hmm_txt) as f:
            text = f.read()

        self.output_gtf = f'{self.outdir}/all_contigs_pfam_a.gtf'
        gtf = open(self.output_gtf, 'w')

        # remove header
        text = text[text.find('\n\n') + 2:]
        # split each query, note here last line is [ok], thus not include
        all_queries = text.split('\n//\n')[:-1]
        # get queries with hits
        queries = []
        for each in all_queries:
            if not '[No hits detected that satisfy reporting thresholds]' in each:
                queries.append(each)

        for query in queries:
            # get the query, accession and description for each query
            line1, line2, line3 = query.split('\n')[:3]
            query_name = line1.split(':')[1].strip()
            query_accession = line2.split(':')[1].strip()
            query_description = line3.split(':')[1].strip()
            # get the contigs between 'Domain annotation for each sequence (and alignments)' and 'Internal pipeline statistics summary'
            text = query.split('Domain annotation for each sequence (and alignments):\n')[1].split(
                'Internal pipeline statistics summary:')[0]
            # get each contigs
            contigs = text.split('>> ')[1:]

            for contig in contigs:
                # extreme case no hit table
                if '[No individual domains that satisfy reporting thresholds (although complete target did)]' in contig:
                    continue
                # get the hit table

                table = contig.split('\n\n  Alignments for each domain:')[0]  # note two space before Alignment!
                line1 = table.split('\n')[0].strip()
                contig_id, frame = line1.split(';frame=')
                frame = int(frame)
                length = int(contig_id.split(';length=')[1])
                # get hits
                hits = table.split('\n')[3:]

                for hit in hits:
                    number, satisfied, score, bias, c_Evalue, i_Evalue, \
                    hmm_from, hmm_to, symbol_1, ali_from, ali_to, symbol_2, env_from, env_to, symbol_3, acc = hit.strip().split()
                    # get '!' hit: satisfy both per-sequence (contig) and per-domain (query) inclusion thresholds
                    if satisfied == '!':
                        start_aa = int(ali_from)
                        end_aa = int(ali_to)
                        i_Evalue = float(i_Evalue)
                        # ali_from and ali_to are aa number, 1-based
                        start_bp = (3 * (start_aa - 1) + 1) + (abs(frame) - 1)
                        end_bp = 3 * end_aa + (abs(frame) - 1)
                        # strand as + or -
                        if frame > 0:
                            strand = '+'
                        else:  # if - strand, ali_from and ali_to counts from the 3', thus need to inverted
                            strand = '-'
                            start_bp, end_bp = length - end_bp + 1, length - start_bp + 1
                        # get attribute
                        attribute = f'name "{query_name}";accession "{query_accession}";description "{query_description}";E_value "{i_Evalue}"'
                        # write to gtf file
                        print(f'{contig_id}\t.\tCDS\t{start_bp}\t{end_bp}\t.\t{strand}\t0\t{attribute}', file=gtf)

        gtf.close()
