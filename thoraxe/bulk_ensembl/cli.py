"""
cli: the ensembl_bulk program (download and index the bulk files).
"""

import argparse
import sys

from thoraxe import utils
from thoraxe.bulk_ensembl import compara, download, index
from thoraxe.bulk_ensembl.layout import BulkDir, DEFAULT_QUERY_SPECIES
from thoraxe.version import __version__


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse ensembl_bulk's command line arguments and returns
    the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog='ensembl_bulk',
        description="""
        ensembl_bulk downloads and indexes the Ensembl bulk files that
        transcript_query reads with --bulk-dir, instead of querying the
        Ensembl REST API and BioMart gene by gene.
        Run "ensembl_bulk download", then "ensembl_bulk index".
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version',
                        version=f'ThorAxe version {__version__}')
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    get = subparsers.add_parser(
        'download', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='download the bulk files of an Ensembl release')
    get.add_argument('bulk_dir', help='bulk folder (created if needed)')
    get.add_argument('-l', '--specieslist', default='',
                     help='restrict the download to these species, as in '
                     'transcript_query: a list separated by commas without '
                     'spaces, or a file with one species per line. If empty, '
                     'every species with an ortholog of the query species is '
                     'downloaded. The query species is always added.')
    get.add_argument('-r', '--release', type=int, default=None,
                     help='Ensembl release, e.g. 115 (required the first time)')
    get.add_argument('-s', '--query-species', default=None,
                     help='species of the genes given to transcript_query '
                     f'(default: {DEFAULT_QUERY_SPECIES})')
    get.add_argument('--homology-species', nargs='*', default=None,
                     help='species whose Compara homology file is also read '
                     '(default: those found for homo_sapiens in release 115; '
                     '"ensembl_bulk index" tells if others are needed)')
    get.add_argument('--biomart-host', default=None,
                     help='BioMart site serving the release (default: its '
                     'Ensembl archive site)')
    get.add_argument('-w', '--workers', type=int, default=4,
                     help='parallel downloads (at most 4)')

    build = subparsers.add_parser(
        'index', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='build the indexes that transcript_query --bulk-dir reads')
    build.add_argument('bulk_dir', help='bulk folder')
    build.add_argument('-l', '--specieslist', default='',
                       help='species to index; all downloaded species if empty')
    build.add_argument('-w', '--workers', type=int, default=1,
                       help='species indexed in parallel (a few GB of memory each)')
    return parser


def main(argv=None):
    "Main function of ensembl_bulk."
    args = parse_command_line().parse_args(argv)
    if args.command == 'download':
        bulk = BulkDir(args.bulk_dir, release=args.release,
                       query_species=args.query_species)
        species_list = utils.species.get_species_list(args.specieslist)
        homology_species = args.homology_species
        if homology_species is None:
            homology_species = compara.DEFAULT_HOMOLOGY_SPECIES.get(
                bulk.query_species, [])
        download.download(bulk, species_list, homology_species,
                          biomart_host=args.biomart_host, workers=args.workers)
    else:
        bulk = BulkDir(args.bulk_dir)
        species_list = utils.species.get_species_list(args.specieslist)
        index.index(bulk, species_list, workers=args.workers)
    return 0


if __name__ == '__main__':
    sys.exit(main())
