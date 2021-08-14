from os import makedirs
from shutil import rmtree
from .MetaGPA import MetaGPA
from .template import Settings
from .tools import get_temp_path


def main(
        control_fq_1: str,
        control_fq_2: str,
        case_fq_1: str,
        case_fq_2: str,
        pfam_hmm: str,
        min_contig_length: int,
        enrichment_cutoff: float,
        outdir: str,
        threads: int,
        debug: bool):

    Main().main(
        control_fq_1=control_fq_1,
        control_fq_2=control_fq_2,
        case_fq_1=case_fq_1,
        case_fq_2=case_fq_2,
        pfam_hmm=pfam_hmm,
        min_contig_length=min_contig_length,
        enrichment_cutoff=enrichment_cutoff,
        outdir=outdir,
        threads=threads,
        debug=debug)


class Main:

    control_fq_1: str
    control_fq_2: str
    case_fq_1: str
    case_fq_2: str
    pfam_hmm: str
    min_contig_length: int
    enrichment_cutoff: float
    outdir: str
    threads: int
    debug: bool

    settings: Settings

    def main(
            self,
            control_fq_1: str,
            control_fq_2: str,
            case_fq_1: str,
            case_fq_2: str,
            pfam_hmm: str,
            min_contig_length: int,
            enrichment_cutoff: float,
            outdir: str,
            threads: int,
            debug: bool):

        self.control_fq_1 = control_fq_1
        self.control_fq_2 = control_fq_2
        self.case_fq_1 = case_fq_1
        self.case_fq_2 = case_fq_2
        self.pfam_hmm = pfam_hmm
        self.min_contig_length = min_contig_length
        self.enrichment_cutoff = enrichment_cutoff
        self.outdir = outdir
        self.threads = threads
        self.debug = debug

        self.set_settings()
        self.makedirs()
        self.execute()
        self.clean_up()

    def set_settings(self):
        self.settings = Settings(
            workdir=get_temp_path(prefix='workdir'),
            outdir=self.outdir,
            threads=self.threads,
            debug=self.debug,
            mock=False)

    def makedirs(self):
        for d in [self.settings.workdir, self.settings.outdir]:
            makedirs(d, exist_ok=True)

    def execute(self):
        MetaGPA(self.settings).main(
            control_fq1=self.control_fq_1,
            control_fq2=self.control_fq_2,
            case_fq1=self.case_fq_1,
            case_fq2=self.case_fq_2,
            pfama=self.pfam_hmm,
            min_contigs_length=self.min_contig_length,
            enrichment_cut_off=self.enrichment_cutoff,
        )

    def clean_up(self):
        if not self.debug:
            rmtree(self.settings.workdir)
