import argparse
import MetaGPA


__version__ = '1.0.0-beta'


PROG = 'python MetaGPA'

DESCRIPTION = f'''\
Metagenomics Genome and Phenome Association (MetaGPA) Analysis

Version: {__version__}

Authors:
  Yu-Cheng Lin (ylin@nycu.edu.edu)
  Weiwei Yang (wyang@neb.com)'''

REQUIRED = [
    {
        'keys': ['-1', '--control-fq-1'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'read 1 fastq file of the control group, .gz format accepted',
        }
    },
    {
        'keys': ['-2', '--control-fq-2'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'read 2 fastq file fo the control group, .gz format accepted',
        }
    },
    {
        'keys': ['-3', '--case-fq-1'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'read 1 fastq file of the case group, .gz format accepted',
        }
    },
    {
        'keys': ['-4', '--case-fq-2'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'read 2 fastq file of the case group, .gz format accepted',
        }
    },
    {
        'keys': ['-p', '--pfam-hmm'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'hmm file for "hmmsearch" algorithm, .gz format accepted',
        }
    },
]

OPTIONAL = [
    {
        'keys': ['-l', '--min-contig-length'],
        'properties': {
            'type': int,
            'required': False,
            'default': 1000,
            'help': 'mininum assembled contig length (default: %(default)s)',
        }
    },
    {
        'keys': ['-c', '--enrichment-cutoff'],
        'properties': {
            'type': float,
            'required': False,
            'default': 3.0,
            'help': 'enrichment score cutoff for the selection of modified contigs (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'MetaGPA-output',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 32,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __version__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        MetaGPA.main(
            control_fq_1=args.control_fq_1,
            control_fq_2=args.control_fq_2,
            case_fq_1=args.case_fq_1,
            case_fq_2=args.case_fq_2,
            pfam_hmm=args.pfam_hmm,
            min_contig_length=args.min_contig_length,
            enrichment_cutoff=args.enrichment_cutoff,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug
        )


if __name__ == '__main__':
    EntryPoint().main()
