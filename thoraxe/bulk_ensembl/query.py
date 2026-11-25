"""
query: transcript_query --bulk-dir, the <gene>/Ensembl/ folder from bulk files.

`run` follows transcript_query's main step by step and writes the same files
with the same writers; only the data source changes. The functions of
transcript_query that filter orthologs, rename the exon table columns and
write tsl.csv, sequences.fasta and errors.csv are reused as they are.
"""

import os
import sys
import warnings
from datetime import datetime

from thoraxe import utils
from thoraxe.transcript_query import transcript_query as tq
from thoraxe.bulk_ensembl import compara, exons, sequences, tsl
from thoraxe.bulk_ensembl.layout import BulkDir, eg_version

#: transcript_query labels the records of the queried gene with this literal
#: text (it is not formatted); kept so that both modes write the same file.
QUERY_EXONS_LABEL = r"{args.species}:{args.genename}"


def _print_if(condition, *args):
    if condition:
        print(*args)


def _check_species(bulk, species_names):
    "Raise if a species has no bulk data or unbuilt indexes."
    problems = {}
    for species in sorted(set(species_names)):
        missing = bulk.missing_indexes(species)
        if missing:
            problems[species] = missing
    if problems:
        details = '; '.join(f"{species}: {', '.join(missing)}"
                            for species, missing in problems.items())
        raise FileNotFoundError(
            f'missing bulk data in {bulk.root} ({details}). Run "ensembl_bulk '
            'download" and "ensembl_bulk index" for these species.')


def resolve_gene(bulk, species, genename, verbose=False):
    "Return the Ensembl gene id of genename (a gene name or an Ensembl id)."
    if tq.is_esemble_id(genename):
        return genename
    print(f"Searching ID for gene with name {genename} in species {species} ...")
    gene_ids = sequences.genes_named(bulk.gtf_db(species), genename)
    _print_if(verbose, f"Found the following list of ids: {gene_ids}")
    if not gene_ids:
        raise Exception(f"No results for {genename}")  # pylint: disable=broad-exception-raised
    print(f"... using gene id {gene_ids[0]} from now on.")
    return gene_ids[0]


def write_ensembl_version(folder, release):
    "Write ensembl_version.csv as transcript_query.save_ensembl_version does."
    with open(os.path.join(folder, 'ensembl_version.csv'), 'w',
              encoding='utf-8') as stream:
        stream.write('Download_date,Download_time,Ensembl_version\n')
        now = datetime.now()
        stream.write(f'{now.strftime("%d-%m-%Y")},{now.strftime("%H:%M")},'
                     f'{eg_version(release)}\n')


def get_orthologs(bulk, gene_id):
    "Bulk version of transcript_query.get_orthologs."
    dortho = compara.homologies(bulk.homologies_db(), gene_id)
    if not dortho:
        raise Exception(  # pylint: disable=broad-exception-raised
            f"Ensembl does not have annotated orthologs for {gene_id}")
    # same species name rule as transcript_query.get_orthologs
    return [d for d in dortho if len(d['target']['species'].split('_')) >= 2]


def run(args):  # pylint: disable=too-many-locals,too-many-statements
    """
    Write <genename>/Ensembl/ from the bulk folder args.bulk_dir.

    args are transcript_query's parsed arguments.
    """
    bulk = BulkDir(args.bulk_dir)
    indexed_for = compara.index_query_species(bulk.homologies_db())
    if indexed_for != args.species:
        raise ValueError(
            f'the homology index of {bulk.root} is for {indexed_for}, not '
            f'{args.species}: run "ensembl_bulk index --query-species '
            f'{args.species}" on a bulk folder downloaded for it.')

    orthokeep = utils.species.get_species_list(args.specieslist)
    if orthokeep is None:
        orthokeep = bulk.ready_species()
        _print_if(args.verbose,
                  f"Using the {len(orthokeep)} species of {bulk.root}")
    _check_species(bulk, [args.species])

    curgene = resolve_gene(bulk, args.species, args.genename, args.verbose)
    cdirectory = args.genename
    _print_if(args.verbose, f"Results will be saved in directory {cdirectory}")
    query_result_subdir = os.path.join(cdirectory, "Ensembl")
    os.makedirs(query_result_subdir, exist_ok=True)

    write_ensembl_version(query_result_subdir, bulk.release)

    print("Writing the gene tree")
    tree_text = compara.tree_newick(bulk.trees_db(), curgene)
    if tree_text is None:
        warnings.warn(f"No gene tree found for gene id {curgene}")
    else:
        with open(os.path.join(query_result_subdir, "tree.nh"), "w",
                  encoding="utf-8") as treeout:
            treeout.write(tree_text)

    print("Looking for orthologs")
    orthologs = get_orthologs(bulk, curgene)
    _print_if(args.verbose, f"Found a total of {len(orthologs)} homologs")
    orthologs_filtered = tq.filter_ortho(orthologs, orthokeep,
                                         relationship=args.orthology)
    _check_species(bulk, [o['target']['species'] for o in orthologs_filtered])

    print("Getting all the transcripts for TSL file")
    tables = tsl.TablesCache(bulk)
    tsl_cur = tables.rows(args.species, curgene)
    tsl_ortho = [tables.rows(o['target']['species'], o['target']['id'])
                 for o in orthologs_filtered]
    tq.write_tsl_file(query_result_subdir, [tsl_cur] + tsl_ortho)

    print("Getting exons sequences")
    # One species open at a time: the genome indexes of all species do not fit in
    # memory (millions of sequence regions for fragmented assemblies). The orthologs
    # come grouped by species, so each genome is still opened once per gene.
    species_sequences = {}

    def exon_records(species, gene_id, label):
        if species not in species_sequences:
            for opened in species_sequences.values():
                opened.close()
            species_sequences.clear()
            species_sequences[species] = sequences.SpeciesSequences(
                bulk.genome_fa(species), bulk.gtf_db(species))
        return species_sequences[species].exon_records(gene_id, label)

    def exon_table(species, gene_id, header):
        if os.path.isfile(bulk.exons_no_dataset(species)):
            return None
        return exons.exon_table(bulk.exons_db(species), gene_id, header=header)

    ffasta = os.path.join(query_result_subdir, "sequences.fasta")
    fexonstable = os.path.join(query_result_subdir, "exonstable.tsv")
    try:
        with open(ffasta, "w", encoding="utf-8") as fastaout, \
                open(fexonstable, "w", encoding="utf-8") as exonstableout:
            extable = exon_table(args.species, curgene, header=True)
            if extable is None:
                tq._store_errors(query_result_subdir, args.species, curgene)  # pylint: disable=protected-access
                sys.exit(1)
            exonstableout.write(tq._rename(extable))  # pylint: disable=protected-access
            for dseq in exon_records(args.species, curgene, QUERY_EXONS_LABEL):
                tq.dictseq2fasta(dseq, dseq['label'], fastaout)
            for ortholog in orthologs_filtered:
                orthoid = ortholog['target']['id']
                orthospecies = ortholog['target']['species']
                ortho_name = f"{orthospecies}:{orthoid}"
                _print_if(args.verbose, f"Getting exons information for {ortho_name}")
                ortho_exontable = exon_table(orthospecies, orthoid, header=False)
                if ortho_exontable is None:
                    warnings.warn(f'No exon table for {orthoid} in {orthospecies}! ')
                    tq._store_errors(query_result_subdir, orthospecies, orthoid)  # pylint: disable=protected-access
                    continue
                exonstableout.write(ortho_exontable)
                for dseq in exon_records(orthospecies, orthoid, ortho_name):
                    tq.dictseq2fasta(dseq, dseq['label'], fastaout)
    finally:
        for opened in species_sequences.values():
            opened.close()
    print("------------------- transcript_query finished -------------------")
