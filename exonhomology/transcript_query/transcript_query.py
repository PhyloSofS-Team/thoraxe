"""
Created on Fri Apr 14 16:32:51 2017

@author: huguesrichard

Un premier jeu de fonctions pour faire des requêtes directement avec
l'API RESTfull de ENSEMBL afin de récupérer tous les gènes homologues
à un gene donné (par son nom courant)
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter

import requests

# REST parameters
SERVER = "https://rest.ensembl.org"

HJSON = {"Content-Type": "application/json"}
HJSONPOST = {"Content-Type": "application/json", "Accept": "application/json"}

NHTREE = {"Content-Type": "text/x-nh"}
BIOMART_HUMAN = "hsapiens_gene_ensembl"

# Other parameters

VERBOSE = False


###
# Utility functions
def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse transcript_query' command line arguments and
    returns the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="transcript_query",
        description="""
        transcript_query download from Ensembl the transcript and exon
        information needed for exonhomology.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genename', type=str, help='gene name in Ensembl (e.g. MAPK8)')
    parser.add_argument(
        '-s',
        '--species',
        help='species to look for the gene name',
        default='homo_sapiens')
    parser.add_argument(
        '-o',
        '--orthology',
        help='Orthology relationship to use; 1:1, 1:n or m:n',
        default='1:n')
    parser.add_argument(
        '-l',
        '--specieslist',
        help='It could be a list of more than one species separated by commas '
        'and without spaces, e.g. homo_sapiens,mus_musculus, or a single file '
        'with the species list (one species per line). If nothing is '
        'indicated, all the available species are used.',
        default='')
    # TO DO: take care of aliases for species names,symbol always use the
    # binomial names when running the code

    return parser.parse_args()


def _check_species_name(species_name):
    """
    Return True if the species_name has the correct format.

    The expected name is binomial but it can also include the subspecies name
    (trinomial). If the species_name doesn't conform the expected format an
    error is raised.

    >>> _check_species_name('homo_sapiens')
    True
    >>> _check_species_name('colobus_angolensis_palliatus')
    True
    >>> _check_species_name('cricetulus_griseus_chok1gshd')
    True
    """
    result = re.match('^[a-z]+_[a-z]+(_[0-9a-z]+)?$', species_name)
    if result is None:
        raise ValueError(
            'The species name should be the binomial/trinomial name lowercased'
            ' with spaces replaced by underscores e.g. Homo sapiens should be '
            'homo_sapiens. {} do not conform the format.'.format(species_name))
    return True


def _read_species_list(file):
    """Return the species list from the file."""
    species = []
    with open(file, 'r') as stream:
        for line in stream:
            species_name = line.strip()
            if _check_species_name(species_name):
                species.append(species_name)
    return species


def _get_species_list(specieslist):
    """
    Return a list of species names (str).

    >>> _get_species_list('')
    >>> _get_species_list('homo_sapiens,mus_musculus')
    ['homo_sapiens', 'mus_musculus']
    """
    if specieslist == '':
        return None

    fields = specieslist.split(',')

    if len(fields) == 1:
        filename = fields[0]
        if os.path.isfile(filename):
            return _read_species_list(filename)
        raise ValueError('{} is not a file.'.format(filename))

    for species_name in fields:
        _check_species_name(species_name)

    return fields


def lodict2csv(listofdicts, out, fnames=None, header=True):
    """
    Write a dictionary list with csv formatting to the stream out.

    fnames: when provided as a list, it is used as the column selection,
        otherwise all keys occuring at least once are used.
    header: should the header be written
    """
    if fnames is None:
        fnames = set([])
        for dictionary in listofdicts:
            fnames.update(list(dictionary.keys()))
        fnames = sorted(fnames)
    csv_writer = csv.DictWriter(
        out, fieldnames=fnames, restval='NA', extrasaction='ignore')
    if header:
        csv_writer.writeheader()
    csv_writer.writerows(listofdicts)
    return len(listofdicts)


def dictseq2fasta(dseq, geneid, out):
    """Write fasta sequences from the exons."""
    colw = 80
    out.write(">%s %s %s\n" % (geneid, dseq['id'], dseq['desc']))
    exseq = dseq['seq']
    exseq = "\n".join([exseq[i:i + colw] for i in range(0, len(exseq), colw)])
    out.write(exseq + "\n")


def species2ensembldataset(species_name):
    """
    Return the name of the ENSEMBL dataset in biomart.

    >>> species2ensembldataset('homo_sapiens')
    'hsapiens_gene_ensembl'
    >>> species2ensembldataset('cebus_capucinus_imitator')
    'ccapucinus_gene_ensembl'
    """
    _check_species_name(species_name)
    names = species_name.split("_")
    biomartname = names[0][0] + names[1] + "_gene_ensembl"
    return biomartname


# Small biomart function from keithshep


def get_biomart_exons_annot(dataset, geneid, header=True):
    """
    Return all transcript information from the dataset and the ensembl geneid.

    ex for MAPK8 in human:
        r = get_biomart_exons_annot("hsapiens_gene_ensembl", "ENSG00000107643")
    Return a requests object (use r.text to get the text for the file)
    """
    # TO DO: Error control on the call

    # see http://ensembl.org/biomart/martview/ for the web application
    biomart_request_url_template = (
        'http://ensembl.org/biomart/martservice?query='
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query  virtualSchemaName = "default" formatter = "TSV" '
        'header = "{ish}" uniqueRows = "0" count = "" '
        'datasetConfigVersion = "0.6">'
        '<Dataset name = "{data}" interface = "default">'
        '<Filter name = "ensembl_gene_id" value = "{eid}"/>'
        '<Attribute name = "ensembl_gene_id" />'
        '<Attribute name = "ensembl_transcript_id" />'
        '<Attribute name = "transcript_start" />'
        '<Attribute name = "transcript_end" />'
        '<Attribute name = "strand" />'
        '<Attribute name = "ensembl_exon_id" />'
        '<Attribute name = "exon_chrom_start" />'
        '<Attribute name = "exon_chrom_end" />'
        '<Attribute name = "rank" />'
        '<Attribute name = "cdna_coding_start" />'
        '<Attribute name = "cdna_coding_end" />'
        '<Attribute name = "genomic_coding_start" />'
        '<Attribute name = "genomic_coding_end" />'
        '<Attribute name = "cds_start" />'
        '<Attribute name = "cds_end" />'
        '<Attribute name = "phase" />'
        '<Attribute name = "end_phase" />'
        '</Dataset>'
        '</Query>')

    biomart_request_url = biomart_request_url_template.format(
        data=dataset, eid=geneid, ish=int(header))
    req = requests.get(biomart_request_url)
    return req


# TO DO : passer toutes les fonctions avec un conteneur générique sur la
# forme de la partie extension


def generic_ensembl_rest_request(extension, params, header):
    """Perform a generic request."""
    request = requests.get(SERVER + extension, params=params, headers=header)
    if VERBOSE:
        print("Trying url:" + request.url)
    if not request.ok:
        request.raise_for_status()
        sys.exit()
    return request


def get_geneids_from_symbol(species, symbol, **params):
    """
    Return gene ID from symbol.

    From a species and a symbol, return the set of geneids corresponding
    to the gene symbol given.
    Uses the /xrefs/symbol RESTfull command
    exple: get_geneids_from_symbol("human", "MAPK8")
    """
    params.setdefault('object_type', 'gene')
    ext_geneid = '/xrefs/symbol/{0}/{1}?'.format(species, symbol)
    request = generic_ensembl_rest_request(ext_geneid, params, HJSON)
    dnames = request.json()
    res = [dictionary["id"] for dictionary in dnames]
    return res


def get_listoftranscripts(ensgeneid, species, **params):
    """
    Return list of transcripts.

    From an ensembl gene id, gets the list of transcripts overlapping
    this gene.
    """
    # TO DO filter on protein coding
    params.setdefault('feature', 'transcript')
    ext_listtr = '/overlap/id/{0}?'.format(ensgeneid)
    request = generic_ensembl_rest_request(ext_listtr, params, HJSON)
    dnames = request.json()
    for dictionary in dnames:
        dictionary['species'] = species
    return dnames


def write_tsl_file(gene_name, l_of_sptr):
    """Write a TSL file from a list of transcripts."""
    fout = "%s_TSL.csv" % (gene_name)
    csvout = open(fout, "w")

    # One trick to get the good names in the header
    other_names = [
        "Species", "Name", "Transcript ID", "Source", "Experiment Source",
        "Biotype", "Flags"
    ]
    cfieldnames = [
        "species", "external_name", "transcript_id", "source", "logic_name",
        "biotype", "transcript_support_level"
    ]
    dnewheader = dict(x for x in zip(cfieldnames, other_names))

    lodict2csv([dnewheader], csvout, fnames=cfieldnames, header=False)
    for sptr in l_of_sptr:
        lodict2csv(sptr, csvout, fnames=cfieldnames, header=False)
    csvout.close()


def get_listofexons(ensgeneid, **params):
    """
    Return list of exons.

    From an ensembl gene id, gets the list of exons that are composing this
    gene by default restricted to the coding exons.
    """
    params.setdefault('feature', 'exon')
    # could think about biotype = protein_coding here
    ext_listex = '/overlap/id/{0}?'.format(ensgeneid)
    request = generic_ensembl_rest_request(ext_listex, params, HJSON)
    dexons = request.json()
    return dexons


def get_exons_sequences(listensexons):  # , **params):
    """
    Return exon sequences.

    From a list of ensembl exons id, it gets the list of exons with their
    sequences.
    """
    # There cannot be more that 50 elements queried at once
    maxitems = 45
    nexons = len(listensexons)
    niter = 1 + nexons // maxitems
    start, end = 0, maxitems
    list_res = []
    for _ in range(niter):
        exons = listensexons[start:end]
        if exons:
            dexons = {"ids": exons}  # , "type": "cds"}
            ext_exons_seq = '/sequence/id/type=cds'
            request = requests.post(
                SERVER + ext_exons_seq,
                headers=HJSONPOST,
                data=json.dumps(dexons))
            # print(request.url)
            if not request.ok:
                print(("FAILED REQUEST: " + str(dexons)))
                request.raise_for_status()
                sys.exit()
            res = request.json()
            list_res = list_res + res
            start, end = end, min(nexons, end + maxitems)
    return list_res


def get_genetree(ensgeneid, **params):
    """
    Return the gene tree.

    Get the gene tree around the gene geneid as of now, the whole tree
    is returned.
    """
    # params.setdefault('object_type',  'gene')
    params.setdefault('nh_format', 'full')
    params.setdefault('aligned', '1')
    ext_genetree = '/genetree/member/id/{0}?'.format(ensgeneid)
    request = generic_ensembl_rest_request(ext_genetree, params, NHTREE)
    return request.text


def get_orthologs(ensgeneid, **params):
    """Get the orthologs from the gene with id ensgeneid."""
    params.setdefault('object_type', 'gene')
    # TO DO: rajouter des parametres pour séparer les différentes relations
    # d'orthologie:
    # -ortholog_one2one,
    # -ortholog_one2many,
    # -within_species_paralog
    # Aussi avoir une option pour filtrer les espèces qui vont être regardées
    # On a besoin de stocker aussi l'info général sur les orthologues
    # pour plus tard
    ext_orthologs = '/homology/id/{0}'.format(ensgeneid)
    request = generic_ensembl_rest_request(ext_orthologs, params, HJSON)
    res = request.json()
    dortho = res['data'][0]['homologies']
    return dortho


def filter_ortho(dortho, species=None, relationship='1:n'):
    """Filter the dictionary of orthologues according to the list of names."""
    # TO DO: rajouter un système de synonymes sur les espèces pour le filtrage
    relationships = _get_relationships(relationship)
    return [
        value for value in dortho if value['type'] in relationships and (
            True if species is None else value['target']['species'].lower() in
            species)
    ]


def get_transcripts_orthologs(ensgeneid, lorthologs):
    """
    Return transcript list from orthologs.

    Wrapper function to call multiple times get_listoftranscripts, given
    a ensembl geneid and a list of orthologs provided by get_orthologs
    Data structure for each ortholog is
    {dn_ds : float, method_link_type : str,rop in
        source : dict(),  target : dict(), taxonomy_level : str,
        type: Enum(ortholog_one2one,
                   ortholog_one2many,
                   within_species_paralog)}
    The dicts for source and target store information about gene sequence,
    data structure:
    {"align_seq" : str, "perc_pos" : float, "id" : str, "protein_id" : str,
     "perc_id" : float, "cigar_line" : str,
     "taxon_id" : int, "species" : str}
    """
    # First get the list for the source
    source_species = lorthologs[0]['source']['species']
    source_transcripts = get_listoftranscripts(ensgeneid, source_species)

    ortho_transcripts = []
    for ortho in lorthologs:
        orthoid = ortho['target']['id']
        orthospecies = ortho['target']['species']
        # orthotaxon = ortho['target']['taxon_id']
        ortho_transcripts.append(get_listoftranscripts(orthoid, orthospecies))

    return source_transcripts, ortho_transcripts


def _get_relationships(notation):
    """
    Return the list of orthology relationships to keep.

    ENSEMBL homology relationships are defined in:
    https://www.ensembl.org/info/genome/compara/homology_types.html

    Posible values to this function are '1:1', '1:n' and 'm:n'.

    >>> _get_relationships('1:1')
    ['ortholog_one2one']
    >>> _get_relationships('1:n')
    ['ortholog_one2one', 'ortholog_one2many']
    >>> _get_relationships('m:n')
    ['ortholog_one2one', 'ortholog_one2many', 'ortholog_many2many']
    """
    if notation == '1:1':
        return ['ortholog_one2one']
    if notation == '1:n':
        return ['ortholog_one2one', 'ortholog_one2many']
    if notation == 'm:n':
        return ['ortholog_one2one', 'ortholog_one2many', 'ortholog_many2many']
    raise ValueError('Orthology should be 1:1, 1:n or m:n')


# TO DO : Refactor main to avoid pylint statements if possible:
# Too many local variables (42/15) and Too many statements (75/50).
def main():  # pylint: disable=too-many-locals,too-many-statements
    """Main script function to download transcript data from ENSEMBL."""
    # 1- Get gene name and species
    # 2- Match with ensembl gene ID
    # 3- Get the set of orthologous genes as ensembl geneID
    # (option for filtering the set of species)
    # 4a- Get all transcript information for the orthologous genes
    # and store the information on the species next to it
    # 4b- Get the exons annotation and sequence
    # 5-  Get the gene tree for the selected species

    # We fix the list of orthologues to a small list inspired from MAPK8
    # orthokeep = [
    #     # "homo_sapiens", "felix_catus", "gallus_gallus",
    #     # "drosophila_melanogaster", "mus_musculus",
    #     # "caenorhabditis_elegans",
    #     # "xenopus_tropicalis", "danio_rerio", "oryctolagus_cuniculus",
    #     # "pan_troglodytes"
    #     'homo_sapiens',
    #     'mus_musculus',
    #     'macaca_mulatta',
    #     'danio_rerio',
    #     'xenopus_tropicalis',
    #     'caenorhabditis_elegans',
    #     'gallus_gallus',
    #     'rattus_norvegicus',
    #     'bos_taurus',
    #     'monodelphis_domestica',
    #     'ornithorhynchus_anatinus',
    #     'drosophila_melanogaster',
    #     'gorilla_gorilla',
    #     'sus_scrofa'
    # ]

    # 1-
    args = parse_command_line()

    # 2-
    orthokeep = _get_species_list(args.specieslist)

    print("Searching ID for gene with name %s in species %s..." %
          (args.genename, args.species))
    geneids = get_geneids_from_symbol(args.species, args.genename)
    print("Found the following list of ids: %s" % (json.dumps(geneids)))
    if not geneids:
        raise KeyError("No gene found, exiting")
    curgene = geneids[0]
    gene_name = "%s_%s" % (args.genename, curgene)
    cdirectory = gene_name
    tsl_subdir = cdirectory + "/TSL/"
    tex_subdir = cdirectory + "/TablesExons/"
    seqsubdir = cdirectory + "/Sequences/"
    print("Using gene id %s from now on." % (curgene))
    print("Results will be saved in directory %s" % (cdirectory))
    if not os.path.exists(cdirectory):
        os.makedirs(tsl_subdir)
        os.makedirs(tex_subdir)
        os.makedirs(seqsubdir)
    # 3-
    # print "Searching for orthologous sequences (ignoring paralogues for now)"
    print("Writing the gene tree")
    tree_text = get_genetree(curgene)
    treeout = open("%s/%s.tree.nh" % (cdirectory, curgene), "w")
    treeout.write(tree_text)
    treeout.close()
    print("Looking for orthologs :")
    orthologs = get_orthologs(curgene)
    nparalogs = len(
        [x for x in orthologs if x['type'] == "within_species_paralog"])
    print("Found a total of %d orthologs, of which %d paralogs" %
          (len(orthologs), nparalogs))
    # ['taxonomy_level']
    print("Orthologous species:")
    number = 0
    corthologs = Counter(
        [ortholog['target']['species'] for ortholog in orthologs])
    for i, j in corthologs.most_common():
        print("  %-23s: %4d" % (i, j))
        # if nt > 5: break
        number += 1
    ##
    orthologs_filtered = filter_ortho(
        orthologs, orthokeep, relationship=args.orthology)
    # TO DO print : orthokeep can be None
    # print("Filtering on %d species, %d matches" % (len(orthokeep),
    #                                                len(orthologs_filtered)))

    print("Getting all the transcripts for preparing a TSL file")
    tsl_cur, tsl_ortho = get_transcripts_orthologs(curgene, orthologs_filtered)
    tsl_out = "%s/%s" % (tsl_subdir, gene_name)
    write_tsl_file(tsl_out, [tsl_cur] + tsl_ortho)

    print("**** Query species : %s" % (args.species))
    print("Got a total of %d transcripts with biotypes" % (len(tsl_cur)))
    for i, j in Counter([dic['biotype'] for dic in tsl_cur]).most_common():
        print("  %-23s: %4d" % (i, j))
    print("**** Orthologues")
    for tr_o in tsl_ortho:
        print("%-22s: %4d transcripts " % (tr_o[0]['species'], len(tr_o)))

    print("Now getting the exons sequences")
    # TO DO revert to multiple files if it is easier
    ffasta = "%s/%s.fasta" % (seqsubdir, gene_name)
    fexonstable = "%s/%s_exonstable.tsv" % (tex_subdir, gene_name)
    fastaout = open(ffasta, "w")
    exonstableout = open(fexonstable, "w")
    dex = get_listofexons(curgene)
    lexid = list({x['exon_id'] for x in dex})
    print("Getting the sequences files for %s" % (curgene))
    curspec_ens_dataset = species2ensembldataset(args.species)
    exfasta = get_exons_sequences(lexid)
    extable = get_biomart_exons_annot(curspec_ens_dataset, curgene)
    exonstableout.write(extable.text)
    exons_name = "%s:%s" % (args.species, args.genename)
    for dseq in exfasta:
        dictseq2fasta(dseq, exons_name, fastaout)
    for ortholog in orthologs_filtered:
        orthoid = ortholog['target']['id']
        orthospecies = ortholog['target']['species']
        # orthotaxon = ortholog['target']['taxon_id']
        ortho_name = "%s:%s" % (orthospecies, orthoid)
        print("Getting exons information for %s" % (ortho_name))
        ortho_ens_dataset = species2ensembldataset(orthospecies)
        dexortho = get_listofexons(orthoid)
        lexidortho = list({x['exon_id'] for x in dexortho})
        print("  - %d exons" % (len(lexidortho)))
        exorthofasta = get_exons_sequences(lexidortho)
        print("  - %d fasta sequences" % (len(exorthofasta)))
        ortho_exontable = get_biomart_exons_annot(
            ortho_ens_dataset, orthoid, header=False)
        print("  - %d lines in the exon table" %
              (ortho_exontable.text.count("\n") + 1))
        exonstableout.write(ortho_exontable.text)
        for dseq in exorthofasta:
            dictseq2fasta(dseq, ortho_name, fastaout)
    fastaout.close()
    exonstableout.close()


if __name__ == '__main__':
    main()
