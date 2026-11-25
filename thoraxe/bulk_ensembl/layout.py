"""
layout: where each bulk file and index lives inside a bulk folder.

::

    <bulk>/bulk_info.json                     release, query species and species list
    <bulk>/compara/emf/Compara.<R>.protein_default.nh.emf.gz
    <bulk>/compara/tsv/<species>/Compara.<R>.protein_default.homologies.tsv.gz
    <bulk>/compara/index/trees.sqlite         gene trees and their members
    <bulk>/compara/index/homologies.sqlite    homologies of the query species
    <bulk>/mysql/<species>_core_<R>_<N>/      core database dump (some tables)
    <bulk>/tables/<species>/                  transcripts.tsv, genes.tsv (tsl.csv)
    <bulk>/exons/raw/<species>.tsv.gz         BioMart dump of the species dataset
    <bulk>/exons/<species>.sqlite             the dump indexed by gene
    <bulk>/gtf/<species>/<Species>.<assembly>.<R>.gtf.gz
    <bulk>/genome/<species>/<Species>.<assembly>.dna.toplevel.fa.gz
    <bulk>/genome_fa/<species>/               decompressed genome + .fai index
    <bulk>/index/<species>.exons.sqlite       GTF genes and exons

``<R>`` is the Ensembl release. Species built by Ensembl Genomes (from release
116: caenorhabditis_elegans, drosophila_melanogaster, saccharomyces_cerevisiae)
also carry its release number ``<EG>``: ``<species>_core_<EG>_<R>_<N>`` and
``<Species>.<assembly>.<EG>.gtf.gz``.
"""

import glob
import json
import os
import re

INFO_FILE = 'bulk_info.json'
DEFAULT_QUERY_SPECIES = 'homo_sapiens'

#: Ensembl Genomes numbering of an Ensembl release, as REST /info/eg_version
#: returns it (release 115 is EG 62).
EG_VERSION_OFFSET = 53


def eg_version(release):
    """
    Return the Ensembl Genomes version that transcript_query records.

    >>> eg_version(115)
    62
    """
    return int(release) - EG_VERSION_OFFSET


def _one(pattern):
    """Return the only file matching pattern."""
    hits = sorted(glob.glob(pattern))
    if len(hits) != 1:
        raise FileNotFoundError(
            f'expected one file matching {pattern}, found {len(hits)}: {hits}')
    return hits[0]


class BulkDir:
    """Paths of one bulk folder (see the module documentation)."""

    def __init__(self, root, release=None, query_species=None):
        self.root = os.path.abspath(root)
        info = self.read_info()
        self.release = int(release or info.get('release') or self._guess_release())
        self.query_species = (query_species or info.get('query_species') or
                              DEFAULT_QUERY_SPECIES)
        self.species = info.get('species')

    # ------------------------------------------------------------- metadata

    @property
    def info_path(self):
        "Path of bulk_info.json."
        return os.path.join(self.root, INFO_FILE)

    def read_info(self):
        "Return the content of bulk_info.json ({} if it does not exist)."
        if not os.path.isfile(self.info_path):
            return {}
        with open(self.info_path, encoding='utf-8') as stream:
            return json.load(stream)

    def write_info(self, **values):
        "Update bulk_info.json with values."
        info = self.read_info()
        info.update(values)
        os.makedirs(self.root, exist_ok=True)
        tmp = self.info_path + '.part'
        with open(tmp, 'w', encoding='utf-8') as stream:
            json.dump(info, stream, indent=1, sort_keys=True)
            stream.write('\n')
        os.replace(tmp, self.info_path)
        self.species = info.get('species')

    def _guess_release(self):
        "Release from the name of the Compara EMF file."
        names = glob.glob(os.path.join(self.root, 'compara', 'emf',
                                       'Compara.*.protein_default.nh.emf.gz'))
        releases = {re.search(r'Compara\.(\d+)\.', name).group(1) for name in names}
        if len(releases) != 1:
            raise ValueError(
                f'cannot tell the Ensembl release of {self.root}: run '
                '"ensembl_bulk download" first or give the release.')
        return releases.pop()

    # ------------------------------------------------------------- compara

    def compara_emf(self):
        "Compara gene trees (Newick EMF dump)."
        return os.path.join(
            self.root, 'compara', 'emf',
            f'Compara.{self.release}.protein_default.nh.emf.gz')

    def compara_tsv(self, species):
        "Compara homologies of one species."
        return os.path.join(
            self.root, 'compara', 'tsv', species,
            f'Compara.{self.release}.protein_default.homologies.tsv.gz')

    def compara_tsv_species(self):
        "Species with a downloaded Compara homology file."
        pattern = self.compara_tsv('*')
        return sorted(os.path.basename(os.path.dirname(path))
                      for path in glob.glob(pattern))

    def trees_db(self):
        "SQLite index of the gene trees."
        return os.path.join(self.root, 'compara', 'index', 'trees.sqlite')

    def homologies_db(self):
        "SQLite index of the homologies of the query species."
        return os.path.join(self.root, 'compara', 'index', 'homologies.sqlite')

    # ------------------------------------------------------------- core dump

    def mysql_dir(self):
        "Folder with the core database dumps."
        return os.path.join(self.root, 'mysql')

    def core_dir(self, species):
        "Core database dump folder of species."
        pattern = os.path.join(self.mysql_dir(), f'{species}_core_*')
        name = (rf'{re.escape(species)}_core_({eg_version(self.release)}_)?'
                rf'{self.release}_\d+')
        hits = [path for path in sorted(glob.glob(pattern))
                if re.fullmatch(name, os.path.basename(path))]
        if len(hits) != 1:
            raise FileNotFoundError(
                f'expected one core dump of {species} for release '
                f'{self.release} in {self.mysql_dir()}, found {hits}')
        return hits[0]

    def tsl_tables_dir(self, species):
        "Transcript and gene tables of species, written by tsl.build_tables."
        return os.path.join(self.root, 'tables', species)

    # ------------------------------------------------------------- BioMart

    def exons_raw(self, species):
        "BioMart dump of species."
        return os.path.join(self.root, 'exons', 'raw', f'{species}.tsv.gz')

    def exons_no_dataset(self, species):
        "Marker written when BioMart has no gene dataset for species."
        return os.path.join(self.root, 'exons', 'raw', f'{species}.no_dataset')

    def exons_db(self, species):
        "SQLite index of the BioMart dump of species."
        return os.path.join(self.root, 'exons', f'{species}.sqlite')

    # ------------------------------------------------------------- genome and GTF

    def gtf_gz(self, species):
        "Complete GTF of species."
        name = rf'[^/]+\.({self.release}|{eg_version(self.release)})\.gtf\.gz'
        hits = [path for path in sorted(glob.glob(os.path.join(
            self.root, 'gtf', species, '*.gtf.gz'))) if re.fullmatch(name, os.path.basename(path))]
        if len(hits) != 1:
            raise FileNotFoundError(
                f'expected one complete GTF of {species} for release {self.release}, '
                f'found {hits}')
        return hits[0]

    def genome_gz(self, species):
        "Unmasked toplevel genome of species."
        return _one(os.path.join(self.root, 'genome', species,
                                 '*.dna.toplevel.fa.gz'))

    def genome_fa(self, species):
        "Decompressed genome of species (indexed with a .fai file)."
        name = os.path.basename(self.genome_gz(species))[:-len('.gz')]
        return os.path.join(self.root, 'genome_fa', species, name)

    def gtf_db(self, species):
        "SQLite index of the GTF genes and exons of species."
        return os.path.join(self.root, 'index', f'{species}.exons.sqlite')

    # ------------------------------------------------------------- readiness

    def missing_indexes(self, species):
        """
        Return the indexes of species that are not built yet (empty if ready).

        A BioMart species without gene dataset counts as ready: its genes get
        recorded in errors.csv, as transcript_query does.
        """
        missing = []
        if not os.path.isfile(os.path.join(self.tsl_tables_dir(species),
                                           'transcripts.tsv')):
            missing.append('tsl tables')
        if not (os.path.isfile(self.exons_db(species)) or
                os.path.isfile(self.exons_no_dataset(species))):
            missing.append('exon table')
        try:
            if not os.path.isfile(self.genome_fa(species) + '.fai'):
                missing.append('genome')
        except FileNotFoundError:
            missing.append('genome')
        if not os.path.isfile(self.gtf_db(species)):
            missing.append('gtf')
        return missing

    def ready_species(self):
        "Species of the bulk folder whose indexes are all built."
        candidates = self.species
        if candidates is None:
            candidates = sorted(
                os.path.basename(os.path.dirname(path)) for path in glob.glob(
                    os.path.join(self.root, 'tables', '*', 'transcripts.tsv')))
        return [species for species in candidates
                if not self.missing_indexes(species)]
