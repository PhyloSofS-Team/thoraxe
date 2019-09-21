"""
Created on Fri Apr 14 16:32:51 2017

@author: huguesrichard & diegozea

Un premier jeu de fonctions pour faire des requêtes directement avec
l'API RESTfull de ENSEMBL afin de récupérer tous les gènes homologues
à un gene donné (par son nom courant)
"""

import argparse
import csv
import json
import os
import sys
import time
import warnings
from collections import Counter
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from thoraxe import utils

# REST parameters
SERVER = "https://rest.ensembl.org"

HJSON = {"Content-Type": "application/json"}
HJSONPOST = {"Content-Type": "application/json", "Accept": "application/json"}

NHTREE = {"Content-Type": "text/x-nh"}
BIOMART_HUMAN = "hsapiens_gene_ensembl"

# Download Utils


def _requests_retry(
        retries=20,
        backoff_factor=1,
        # github.com/Ensembl/ensembl-rest/wiki/HTTP-Response-Codes
        # OLD: (500, 501, 502, 503, 504, 401, 403, 404),
        status_forcelist=(403, 408, 503),
        session=None):
    """
    Request with retry and wait.

    Code taking from Sandilya Bhamidipati blog post at
    https://dev.to/ssbozy/python-requests-with-retries-4p03
    """
    if session is None:
        session = requests.Session()

    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('https://', adapter)

    return session


SESSION = _requests_retry()


def _request_ensembl_retry(session, session_method, *args, **kargs):
    """
    Try to request Ensembl waiting if needed.

    https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
    """
    response = getattr(session, session_method)(*args, **kargs)
    wait = response.headers.get('Retry-After')
    if wait:
        warnings.warn(
            'Ensembl rate limit reached, waiting for {} seconds.'.format(wait))
        time.sleep(float(wait) + 1.0)
        response = getattr(session, session_method)(*args, **kargs)

    return response


def _check_biomart_response(response):
    """
    Return True if the response doesn't look like an query ERROR or 404 page.

    >>> _check_biomart_response('')
    False
    >>> _check_biomart_response('<html>')
    False
    >>> _check_biomart_response('Query ERROR: caught BioMart::Exception')
    False
    >>> _check_biomart_response('ENSG00000107643	ENST00000432379	48306655')
    True
    """
    response = response.strip()
    if response.startswith('<html') or response.startswith(
            'Query ERROR:') or not response:
        return False

    return True


def _request_ensembl_redirect(*args, **kargs):
    """Try to request Ensembl with and without redirect."""
    response = _request_ensembl_retry(SESSION, "get", *args, **kargs)

    if response.ok:
        return response

    response = _request_ensembl_retry(requests,
                                      "get",
                                      *args,
                                      **kargs,
                                      allow_redirects=False)
    if response.ok:
        return response

    warnings.warn('Failed request for args: {} kargs: {}'.format(args, kargs))
    return None


def generic_ensembl_rest_request(extension, params, header):
    "Perform a generic request."
    return _request_ensembl_redirect(SERVER + extension,
                                     params=params,
                                     headers=header)


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
        information needed for thoraxe.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('genename',
                        type=str,
                        help='gene name in Ensembl (e.g. MAPK8)')
    parser.add_argument('-s',
                        '--species',
                        help='species to look for the gene name',
                        default='homo_sapiens')
    parser.add_argument('-o',
                        '--orthology',
                        help='Orthology relationship to use; 1:1, 1:n or m:n',
                        default='1:1')
    parser.add_argument(
        '-l',
        '--specieslist',
        help='It could be a list of more than one species separated by commas '
        'and without spaces, e.g. homo_sapiens,mus_musculus, or a single file '
        'with the species list (one species per line). If nothing is '
        'indicated, all the available species are used.',
        default='')
    parser.add_argument('-v',
                        '--verbose',
                        help='Print detailed progress.',
                        action='store_true')
    # TO DO: take care of aliases for species names,symbol always use the
    # binomial names when running the code

    return parser.parse_args()


def lodict2csv(listofdicts, out, fnames=None, header=True):
    """
    Write a dictionary list with csv formatting to the stream out.

    :param fnames: when provided as a list, it is used as the column selection,
                   otherwise all keys occuring at least once are used.
    :param header: should the header be written
    """
    if fnames is None:
        fnames = set([])
        for dictionary in listofdicts:
            fnames.update(list(dictionary.keys()))
        fnames = sorted(fnames)
    csv_writer = csv.DictWriter(out,
                                fieldnames=fnames,
                                restval='NA',
                                extrasaction='ignore')
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


def _species2ensembldataset(species_name):
    """
    Return the name of the ENSEMBL dataset in biomart.

    >>> _species2ensembldataset('homo_sapiens')
    ['hsapiens_gene_ensembl']
    >>> _species2ensembldataset('cebus_capucinus_imitator')
    ['ccapucinus_gene_ensembl', 'cimitator_gene_ensembl']
    """
    utils.species.check_species_name(species_name)
    names = species_name.split("_")
    return [
        names[0][0] + names[i] + "_gene_ensembl" for i in range(1, len(names))
    ]


# Small biomart function from keithshep


def _biomart_exons_annot_request(dataset, geneid, header=True):
    """
    Return all transcript information from the dataset and the ensembl geneid.

    For example, MAPK8 in human:
    _biomart_exons_annot_request("hsapiens_gene_ensembl", "ENSG00000107643")
    """
    # TO DO: Error control on the call

    # see http://ensembl.org/biomart/martview/ for the web application
    biomart_request_url_template = (
        'https://ensembl.org/biomart/martservice?query='
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query  virtualSchemaName = "default" formatter = "TSV" '
        'header = "{ish}" uniqueRows = "0" count = "" '
        'datasetConfigVersion = "0.6">'
        '<Dataset name = "{data}" interface = "default">'
        '<Filter name = "ensembl_gene_id" value = "{eid}"/>'
        '<Attribute name = "ensembl_gene_id" />'
        '<Attribute name = "ensembl_transcript_id" />'
        '<Attribute name = "ensembl_peptide_id" />'
        '<Attribute name = "strand" />'
        '<Attribute name = "ensembl_exon_id" />'
        '<Attribute name = "exon_chrom_start" />'
        '<Attribute name = "exon_chrom_end" />'
        '<Attribute name = "rank" />'
        '<Attribute name = "cdna_coding_start" />'
        '<Attribute name = "cdna_coding_end" />'
        '<Attribute name = "genomic_coding_start" />'
        '<Attribute name = "genomic_coding_end" />'
        '<Attribute name = "phase" />'
        '<Attribute name = "end_phase" />'
        '</Dataset>'
        '</Query>')

    biomart_request_url = biomart_request_url_template.format(data=dataset,
                                                              eid=geneid,
                                                              ish=int(header))
    try:
        req = _request_ensembl_redirect(biomart_request_url)
        if req.status_code >= 300:
            warnings.warn('BioMart request status for {} in {}: {}.'.format(
                geneid, dataset, req.status_code))
        response = req.text
    except Exception as err:  # pylint: disable=broad-except
        print(err)
        print('Query:\n', biomart_request_url, '\n')
        response = ''

    return response


def get_biomart_exons_annot(species_name, geneid, header=True):
    """Return transcript information from a ensembl geneid and species name."""
    dataset_names = _species2ensembldataset(species_name)
    # Try the download option quickly before retrying:
    for dataset in dataset_names:
        response = _biomart_exons_annot_request(dataset, geneid, header=header)
        if _check_biomart_response(response):
            return response

    raise Exception(
        'It can not found {} in biomart (tried: {}).\nLast response:\n{}'.
        format(species_name, dataset_names, response))


# TO DO : passer toutes les fonctions avec un conteneur générique sur la
# forme de la partie extension


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


def write_tsl_file(path, l_of_sptr):
    """Write a TSL file from a list of transcripts."""
    with open(os.path.join(path, "tsl.csv"), "w") as csvout:
        # One trick to get the good names in the header
        other_names = [
            "Species", "Name", "TranscriptID", "Source", "ExperimentSource",
            "Biotype", "Flags", "Version"
        ]
        cfieldnames = [
            "species", "external_name", "transcript_id", "source",
            "logic_name", "biotype", "transcript_support_level", "version"
        ]
        dnewheader = dict(x for x in zip(cfieldnames, other_names))

        lodict2csv([dnewheader], csvout, fnames=cfieldnames, header=False)

        for sptr in l_of_sptr:
            lodict2csv(sptr, csvout, fnames=cfieldnames, header=False)


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
            request = _request_ensembl_retry(SESSION,
                                             "post",
                                             SERVER + ext_exons_seq,
                                             headers=HJSONPOST,
                                             data=json.dumps(dexons))

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
    params.setdefault('object_type', 'gene')
    params.setdefault('nh_format', 'full')  # 'species'
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


FORBIDDEN_SPECIES = {'bos_taurus_hybrid', 'bos_indicus_hybrid'}


def _allowed_species(species):
    "True if the species name is not in FORBIDDEN_SPECIES"
    return species not in FORBIDDEN_SPECIES


def filter_ortho(dortho, species=None, relationship='1:n'):
    """Filter the dictionary of orthologues according to the list of names."""
    # TO DO: rajouter un système de synonymes sur les espèces pour le filtrage
    relationships = _get_relationships(relationship)
    orthologs = [
        value for value in dortho if value['type'] in relationships and (
            _allowed_species(value['target']['species'].lower()) and
            (True if species is None else value['target']['species'].lower() in
             species))
    ]
    if dortho and not orthologs:
        raise Exception(
            "There are not {} in the requested species.".format(relationships))
    return orthologs


def get_transcripts_orthologs(ensgeneid, lorthologs):
    """
    Return transcript list from orthologs.

    Wrapper function to call multiple times get_listoftranscripts, given
    a ensembl geneid and a list of orthologs provided by get_orthologs
    Data structure for each ortholog is
    ::

        {dn_ds : float, method_link_type : str,rop in
            source : dict(),  target : dict(), taxonomy_level : str,
            type: Enum(ortholog_one2one,
                       ortholog_one2many,
                       within_species_paralog)}

    The dicts for source and target store information about gene sequence,
    data structure:
    ::

        {"align_seq" : str, "perc_pos" : float, "id" : str,
        "protein_id" : str, "perc_id" : float, "cigar_line" : str,
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


def _rename(table_text):
    """Rename columns"""
    table_text = table_text.replace('Gene stable ID', 'GeneID', 1)
    table_text = table_text.replace('Transcript stable ID', 'TranscriptID', 1)
    table_text = table_text.replace('Protein stable ID', 'ProteinID', 1)
    table_text = table_text.replace('Exon stable ID', 'ExonID', 1)
    table_text = table_text.replace('Exon region start (bp)',
                                    'ExonRegionStart', 1)
    table_text = table_text.replace('Exon region end (bp)', 'ExonRegionEnd', 1)
    table_text = table_text.replace('Exon rank in transcript', 'ExonRank', 1)
    table_text = table_text.replace('cDNA coding start', 'cDNA_CodingStart', 1)
    table_text = table_text.replace('cDNA coding end', 'cDNA_CodingEnd', 1)
    table_text = table_text.replace('Genomic coding start',
                                    'GenomicCodingStart', 1)
    table_text = table_text.replace('Genomic coding end', 'GenomicCodingEnd',
                                    1)
    table_text = table_text.replace('Start phase', 'StartPhase', 1)
    table_text = table_text.replace('End phase', 'EndPhase', 1)
    return table_text


def _print_if(condition, *args):
    """Print only if condition is True."""
    if condition:
        print(*args)


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

    # 1-
    args = parse_command_line()

    # 2-
    orthokeep = utils.species.get_species_list(args.specieslist)

    print("Searching ID for gene with name %s in species %s ..." %
          (args.genename, args.species))
    geneids = get_geneids_from_symbol(args.species, args.genename)
    _print_if(args.verbose,
              "Found the following list of ids: %s" % (json.dumps(geneids)))
    if not geneids:
        raise Exception("No results for {}".format(args.genename))
    curgene = geneids[0]
    gene_name = args.genename
    cdirectory = gene_name
    query_result_subdir = os.path.join(cdirectory, "Ensembl")
    print("... using gene id %s from now on." % (curgene))
    _print_if(args.verbose,
              "Results will be saved in directory %s" % (cdirectory))
    if not os.path.exists(cdirectory):
        os.makedirs(query_result_subdir)

    # 3-
    # print "Searching for orthologous sequences (ignoring paralogues for now)"
    print("Writing the gene tree")
    tree_text = get_genetree(curgene)
    with open(os.path.join(query_result_subdir, "tree.nh"), "w") as treeout:
        treeout.write(tree_text)

    print("Looking for orthologs")
    orthologs = get_orthologs(curgene)
    nparalogs = len(
        [x for x in orthologs if x['type'] == "within_species_paralog"])
    _print_if(
        args.verbose, "Found a total of %d orthologs, of which %d paralogs" %
        (len(orthologs), nparalogs))
    # ['taxonomy_level']
    _print_if(args.verbose, "Orthologous species:")
    number = 0
    corthologs = Counter(
        [ortholog['target']['species'] for ortholog in orthologs])
    for i, j in corthologs.most_common():
        _print_if(args.verbose, "  %-23s: %4d" % (i, j))
        # if nt > 5: break
        number += 1
    ##
    orthologs_filtered = filter_ortho(orthologs,
                                      orthokeep,
                                      relationship=args.orthology)
    # TO DO print : orthokeep can be None
    # _print_if(args.verbose,
    # "Filtering on %d species, %d matches" % (len(orthokeep),
    #                                          len(orthologs_filtered)))

    print("Getting all the transcripts for TSL file")
    tsl_cur, tsl_ortho = get_transcripts_orthologs(curgene, orthologs_filtered)
    write_tsl_file(query_result_subdir, [tsl_cur] + tsl_ortho)

    _print_if(args.verbose, "**** Query species : %s" % (args.species))
    _print_if(args.verbose,
              "Got a total of %d transcripts with biotypes" % (len(tsl_cur)))
    for i, j in Counter([dic['biotype'] for dic in tsl_cur]).most_common():
        _print_if(args.verbose, "  %-23s: %4d" % (i, j))
    _print_if(args.verbose, "**** Orthologues")
    for tr_o in tsl_ortho:
        _print_if(args.verbose,
                  "%-22s: %4d transcripts " % (tr_o[0]['species'], len(tr_o)))

    print("Getting exons sequences")
    # TO DO revert to multiple files if it is easier
    ffasta = os.path.join(query_result_subdir, "sequences.fasta")
    fexonstable = os.path.join(query_result_subdir, "exonstable.tsv")
    fastaout = open(ffasta, "w")
    exonstableout = open(fexonstable, "w")
    dex = get_listofexons(curgene)
    lexid = list({x['exon_id'] for x in dex})
    _print_if(args.verbose, "Getting the sequences files for %s" % (curgene))
    exfasta = get_exons_sequences(lexid)
    extable = get_biomart_exons_annot(args.species, curgene)
    extable = _rename(extable)
    exonstableout.write(extable)
    exons_name = "%s:%s" % (args.species, args.genename)
    for dseq in exfasta:
        dictseq2fasta(dseq, exons_name, fastaout)
    for ortholog in orthologs_filtered:
        orthoid = ortholog['target']['id']
        orthospecies = ortholog['target']['species']
        # orthotaxon = ortholog['target']['taxon_id']
        ortho_name = "%s:%s" % (orthospecies, orthoid)
        _print_if(args.verbose,
                  "Getting exons information for %s" % (ortho_name))
        dexortho = get_listofexons(orthoid)
        lexidortho = list({x['exon_id'] for x in dexortho})
        _print_if(args.verbose, "  - %d exons" % (len(lexidortho)))
        exorthofasta = get_exons_sequences(lexidortho)
        _print_if(args.verbose, "  - %d fasta sequences" % (len(exorthofasta)))
        ortho_exontable = get_biomart_exons_annot(orthospecies,
                                                  orthoid,
                                                  header=False)
        _print_if(
            args.verbose, "  - %d lines in the exon table" %
            (ortho_exontable.count("\n") + 1))
        exonstableout.write(ortho_exontable)
        for dseq in exorthofasta:
            dictseq2fasta(dseq, ortho_name, fastaout)
    fastaout.close()
    exonstableout.close()


if __name__ == '__main__':
    main()
