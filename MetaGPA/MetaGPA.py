from .mapping import Mapping
from .annotation import Annotation
from .selection import SelectContigs
from .plot import PlotContigEnrichment
from .template import Processor, Settings
from .enrichment import CalculateEnrichment
from .statistics import CorrectedFisherTest
from .assembly import Assembly, CombineAssembly
from .count_pfam import CountNonredundantPfams, MergePfamCounts


class MetaGPA(Processor):

    control_fq1: str
    control_fq2: str
    case_fq1: str
    case_fq2: str
    pfama: str
    min_contigs_length: int
    enrichment_cut_off: float

    control_fa: str
    case_fa: str
    all_fa: str
    all_nd_fa: str
    gtf: str
    control_bam: str
    case_bam: str
    enrichment_csv: str
    modified_gtf: str
    unmodified_gtf: str
    modified_pfam_count_csv: str
    unmodified_pfam_count_csv: str
    merged_pfam_count_csv: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            control_fq1: str,
            control_fq2: str,
            case_fq1: str,
            case_fq2: str,
            pfama: str,
            min_contigs_length: int,
            enrichment_cut_off: float):

        self.control_fq1 = control_fq1
        self.control_fq2 = control_fq2
        self.case_fq1 = case_fq1
        self.case_fq2 = case_fq2
        self.pfama = pfama
        self.min_contigs_length = min_contigs_length
        self.enrichment_cut_off = enrichment_cut_off

        self.assembly()
        self.combine_assembly()
        self.remove_redundant_contigs()
        self.annotation()
        self.mapping()
        self.calculate_enrichment()
        self.plot_contig_enrichment()
        self.select_contigs()
        self.count_nonredundant_pfams()
        self.merge_pfam_counts()
        self.corrected_fisher_test()

    def assembly(self):
        assembly = Assembly(self.settings).main

        self.control_fa = assembly(
            fq1=self.control_fq1,
            fq2=self.control_fq2,
            assembly_name='control',
            min_contigs_length=self.min_contigs_length)

        self.case_fa = assembly(
            fq1=self.case_fq1,
            fq2=self.case_fq2,
            assembly_name='case',
            min_contigs_length=self.min_contigs_length)

    def combine_assembly(self):
        self.all_fa = CombineAssembly(self.settings).main(
            control_fa=self.control_fa,
            case_fa=self.case_fa)

    def remove_redundant_contigs(self):
        self.all_nd_fa = f'{self.outdir}/all_contigs.nd.fa'
        cmd = f'cd-hit-est -i {self.all_fa} -o {self.all_nd_fa} -c 0.95 -n 10 -d 0 -M 0 -T {self.threads}'
        self.call(cmd)

    def annotation(self):
        self.gtf = Annotation(self.settings).main(
            fna=self.all_nd_fa,
            pfama=self.pfama)

    def mapping(self):
        self.control_bam, self.case_bam = Mapping(self.settings).main(
            fna=self.all_nd_fa,
            control_fq1=self.control_fq1,
            control_fq2=self.control_fq2,
            case_fq1=self.case_fq1,
            case_fq2=self.case_fq2)

    def calculate_enrichment(self):
        self.enrichment_csv = CalculateEnrichment(self.settings).main(
            fna=self.all_nd_fa,
            control_fq1=self.control_fq1,
            case_fq1=self.case_fq1,
            control_bam=self.control_bam,
            case_bam=self.case_bam)

    def plot_contig_enrichment(self):
        PlotContigEnrichment(self.settings).main(
            enrichment_csv=self.enrichment_csv)

    def select_contigs(self):
        self.modified_gtf, self.unmodified_gtf = SelectContigs(self.settings).main(
            fna=self.all_nd_fa,
            gtf=self.gtf,
            enrichment_csv=self.enrichment_csv,
            enrichment_cut_off=self.enrichment_cut_off)

    def count_nonredundant_pfams(self):
        self.modified_pfam_count_csv = f'{self.outdir}/modified_contigs_pfam_a_counts.csv'
        self.unmodified_pfam_count_csv = f'{self.outdir}/unmodified_contigs_pfam_a_counts.csv'

        for gtf, csv in [
            (self.modified_gtf, self.modified_pfam_count_csv),
            (self.unmodified_gtf, self.unmodified_pfam_count_csv)
        ]:
            CountNonredundantPfams(self.settings).main(
                input_gtf=gtf, output_csv=csv)

    def merge_pfam_counts(self):
        self.merged_pfam_count_csv = MergePfamCounts(self.settings).main(
            modified_pfam_count_csv=self.modified_pfam_count_csv,
            unmodified_pfam_count_csv=self.unmodified_pfam_count_csv)

    def corrected_fisher_test(self):
        CorrectedFisherTest(self.settings).main(
            merged_pfam_count_csv=self.merged_pfam_count_csv)
