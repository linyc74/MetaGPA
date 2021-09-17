import os
import shutil
import unittest
from MetaGPA.MetaGPA import MetaGPA
from MetaGPA.template import Settings


class TestMetaGPA(unittest.TestCase):

    def setUp(self):
        self.set_dirs()
        self.makedirs()
        self.set_settings()

    def set_dirs(self):
        self.basedir = os.path.dirname(__file__)
        self.workdir = f'{self.basedir}/workdir'
        self.outdir = f'{self.basedir}/outdir'

    def makedirs(self):
        for d in [self.workdir, self.outdir]:
            os.makedirs(d, exist_ok=True)

    def set_settings(self):
        self.settings = Settings(
            workdir=self.workdir,
            outdir=self.outdir,
            threads=4,
            memory=8,
            debug=True,
            mock=False)

    def tearDown(self):
        for d in [self.workdir, self.outdir]:
            shutil.rmtree(d)

    def test_main(self):
        MetaGPA(self.settings).main(
            control_fq1=f'{self.basedir}/test_control.1.fq.gz',
            control_fq2=f'{self.basedir}/test_control.2.fq.gz',
            case_fq1=f'{self.basedir}/test_case.1.fq.gz',
            case_fq2=f'{self.basedir}/test_case.2.fq.gz',
            pfama=f'{self.basedir}/Pfam-A-sampled.hmm.gz',
            min_contigs_length=1,
            enrichment_cut_off=1.0
        )
