#!/usr/bin/env python3
"""
Created on Fri Apr 14 16:32:51 2017

@author: huguesrichard

Un premier jeu de fonctions pour faire des requêtes directement avec
l'API RESTfull de ENSEMBL afin de récupérer tous les gènes homologues
à un gene donné (par son nom courant)
"""

import csv
import io
import json
import os
import sys
from collections import Counter

import requests

# REST parameters
server = "https://rest.ensembl.org"

hjson = {"Content-Type": "application/json"}
hjsonpost = {"Content-Type": "application/json", "Accept": "application/json"}

nhtree = {"Content-Type": "text/x-nh"}
biomartHuman = "hsapiens_gene_ensembl"

# Other parameters

verbose = False


###
# Utility functions
def lodict2csvstring(listofdicts):
    si = io.BytesIO()
    fnames = set([])
    for d in listofdicts:
        fnames.update(list(d.keys()))
    fnames = sorted(fnames)
    cw = csv.DictWriter(si, fieldnames=fnames, restval='NA')
    cw.writeheader()
    cw.writerows(listofdicts)
    return si.getvalue()


def lodict2csv(listofdicts, out, fnames=None, header=True):
    """
    Write a list of dictionnary with csv formatting to the stream out.

    fnames: when provided as a list, it is used as the column selection,
        otherwise all keys occuring at least once are used.
    header: should the header be written
    """
    if fnames is None:
        fnames = set([])
        for d in listofdicts:
            fnames.update(list(d.keys()))
        fnames = sorted(fnames)
    cw = csv.DictWriter(
        out, fieldnames=fnames, restval='NA', extrasaction='ignore')
    if header:
        cw.writeheader()
    cw.writerows(listofdicts)
    return len(listofdicts)


def dictseq2fasta(dseq, geneid, out):
    """Write fasta sequences from the exons."""
    colw = 80
    out.write(">%s %s %s\n" % (geneid, dseq['id'], dseq['desc']))
    exseq = dseq['seq']
    exseq = "\n".join([exseq[i:i + colw] for i in range(0, len(exseq), colw)])
    out.write(exseq + "\n")


def species2ensembldataset(speciesname):
    binname = speciesname.split("_")
    assert len(binname) == 2, "Error for the species name"
    biomartname = binname[0][0] + binname[1] + "_gene_ensembl"
    return biomartname


# Small biomart function from keithshep


def get_biomart_exons_annot(dataset, geneid, header=True):
    """
    Return all transcript information from the dataset and the ensembl geneid.

    ex for MAPK8 in human:
        r = get_biomart_exons_annot("hsapiens_gene_ensembl", "ENSG00000107643")
    Return a requests object (use r.text to get the text for the file)
    """
    # TODO: Error control on the call
    ##
    # see http://ensembl.org/biomart/martview/ for the web application
    biomart_request_url_template = """
    http://ensembl.org/biomart/martservice?query=
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "{ish}"
        uniqueRows = "0" count = "" datasetConfigVersion = "0.6">
    <Dataset name = "{data}" interface = "default">
    <Filter name = "ensembl_gene_id" value = "{eid}"/>
    <Attribute name = "ensembl_gene_id" />
    <Attribute name = "ensembl_transcript_id" />
    <Attribute name = "transcript_start" />
    <Attribute name = "transcript_end" />
    <Attribute name = "strand" />
    <Attribute name = "ensembl_exon_id" />
    <Attribute name = "exon_chrom_start" />
    <Attribute name = "exon_chrom_end" />
    <Attribute name = "rank" />
    <Attribute name = "cdna_coding_start" />
    <Attribute name = "cdna_coding_end" />
    <Attribute name = "genomic_coding_start" />
    <Attribute name = "genomic_coding_end" />
    <Attribute name = "cds_start" />
    <Attribute name = "cds_end" />
    <Attribute name = "phase" />
    <Attribute name = "end_phase" />
    </Dataset>
    </Query>
    """

    biomart_request_url = biomart_request_url_template.format(
        data=dataset, eid=geneid, ish=int(header))
    req = requests.get(biomart_request_url)
    return req


# TODO : passer toutes les fonctions avec un conteneur générique sur la
# forme de la partie extension


def generic_ensemblREST_request(extension, params, header):
    """Perform a generic request."""
    r = requests.get(server + extension, params=params, headers=header)
    if verbose:
        print("Trying url:" + r.url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r


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
    r = generic_ensemblREST_request(ext_geneid, params, hjson)
    dnames = r.json()
    res = [x["id"] for x in dnames]
    return res


def get_listoftranscripts(ensgeneid, species, **params):
    """
    Return list of transcripts.

    From an ensembl gene id, gets the list of transcripts overlapping
    this gene.
    """
    # TODO filter on protein coding
    params.setdefault('feature', 'transcript')
    ext_listtr = '/overlap/id/{0}?'.format(ensgeneid)
    r = generic_ensemblREST_request(ext_listtr, params, hjson)
    dnames = r.json()
    for x in dnames:
        x['species'] = species
    return (dnames)


def write_TSL_file(gene_name, l_of_sptr):
    """Write a TSL file from a list of transcripts."""
    fout = "%s_TSL.csv" % (gene_name)
    csvout = open(fout, "w")

    # One trick to get the good names in the header
    OtherNames = [
        "Species", "Name", "Transcript ID", "Source", "Experiment Source",
        "Biotype", "Flags"
    ]
    cfieldnames = [
        "species", "external_name", "transcript_id", "source", "logic_name",
        "biotype", "transcript_support_level"
    ]
    dnewheader = dict(x for x in zip(cfieldnames, OtherNames))

    lodict2csv([dnewheader], csvout, fnames=cfieldnames, header=False)
    for sptr in l_of_sptr:
        lodict2csv(sptr, csvout, fnames=cfieldnames, header=False)
    csvout.close()
    return None


def get_listofexons(ensgeneid, **params):
    """
    Return list of exons.

    From an ensembl gene id, gets the list of exons that are composing this
    gene by default restricted to the coding exons.
    """
    params.setdefault('feature', 'exon')
    # could think about biotype = protein_coding here
    ext_listex = '/overlap/id/{0}?'.format(ensgeneid)
    r = generic_ensemblREST_request(ext_listex, params, hjson)
    dexons = r.json()
    return dexons


def get_exons_sequences(listensexons, **params):
    """
    Return exon sequences.

    From a list of ensembl exons id, it gets the list of exons with their
    sequences.
    """
    # There cannot be more that 50 elements queried at once
    MAXITEMS = 45
    nexons = len(listensexons)
    niter = 1 + nexons // MAXITEMS
    s, e = 0, MAXITEMS
    list_res = []
    for i in range(niter):
        exons = listensexons[s:e]
        if len(exons) > 0:
            dexons = {"ids": exons}  # , "type": "cds"}
            ext_exons_seq = '/sequence/id/type=cds'
            r = requests.post(
                server + ext_exons_seq,
                headers=hjsonpost,
                data=json.dumps(dexons))
            # print(r.url)
            if not r.ok:
                print(("FAILED REQUEST: " + str(dexons)))
                r.raise_for_status()
                sys.exit()
            res = r.json()
            list_res = list_res + res
            s, e = e, min(nexons, e + MAXITEMS)
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
    r = generic_ensemblREST_request(ext_genetree, params, nhtree)
    return r.text


def get_orthologs(ensgeneid, **params):
    """Get the orthologs from the gene with id ensgeneid."""
    params.setdefault('object_type', 'gene')
    # TODO: rajouter des parametres pour séparer les différentes relations
    # d'orthologie:
    # -ortholog_one2one,
    # -ortholog_one2many,
    # -within_species_paralog
    # Aussi avoir une option pour filtrer les espèces qui vont être regardées
    # On a besoin de stocker aussi l'info général sur les orthologues
    # pour plus tard
    ext_orthologs = '/homology/id/{0}'.format(ensgeneid)
    r = generic_ensemblREST_request(ext_orthologs, params, hjson)
    res = r.json()
    dortho = res['data'][0]['homologies']
    return dortho


def filter_ortho(l_ortho,
                 dortho,
                 relation=[
                     "ortholog_one2one", "ortholog_one2many",
                     "within_species_paralog"
                 ]):
    """Filter the dictionary of orthologues according to the list of names."""
    # TODO: rajouter un système de synonymes sur les espèces pour
    # le filtrage
    filtered_list = [
        v for v in dortho
        if v['target']['species'].lower() in l_ortho and v['type'] in relation
    ]
    return filtered_list


def get_transcripts_orthologs(ensgeneid, lorthologs):
    """
    Return transcript list from orthologs.

    Wrapper function to call multiple times get_listoftranscripts, given
    a ensembl geneid and a list of orthologs provided by get_orthologs
    Data structure for each ortholog is
    {dn_ds : float, method_link_type : str,
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
    for o in lorthologs:
        orthoid = o['target']['id']
        orthospecies = o['target']['species']
        # orthotaxon = o['target']['taxon_id']
        ortho_transcripts.append(get_listoftranscripts(orthoid, orthospecies))

    return source_transcripts, ortho_transcripts


def main():
    # 1- Get gene name and species
    # 2- Match with ensembl gene ID
    # 3- Get the set of orthologous genes as ensembl geneID
    # (option for filtering the set of species)
    # 4a- Get all transcript information for the orthologous genes
    # and store the information on the species next to it
    # 4b- Get the exons annotation and sequence
    # 5-  Get the gene tree for the selected species

    # We fix the list of orthologues to a small list inspired from MAPK8
    orthokeep = [
        # "homo_sapiens", "felix_catus", "gallus_gallus",
        # "drosophila_melanogaster", "mus_musculus", "caenorhabditis_elegans",
        # "xenopus_tropicalis", "danio_rerio", "oryctolagus_cuniculus",
        # "pan_troglodytes"
        'homo_sapiens',
        'mus_musculus',
        'macaca_mulatta',
        'danio_rerio',
        'xenopus_tropicalis',
        'caenorhabditis_elegans',
        'gallus_gallus',
        'rattus_norvegicus',
        'bos_taurus',
        'monodelphis_domestica',
        'ornithorhynchus_anatinus',
        'drosophila_melanogaster',
        'gorilla_gorilla',
        'sus_scrofa'
    ]

    # 1-
    if len(sys.argv) == 3:
        species, symbol = sys.argv[1:]
    elif len(sys.argv) == 2:
        species = 'homo_sapiens'
        symbol = sys.argv[1]
    else:
        print("No parameter given, running with default")
        # species, symbol = 'homo_sapiens', 'DNM2'
        species, symbol = 'homo_sapiens', 'MAPK8'
        # species, symbol = 'homo_sapiens', 'snap25'
        # species, symbol = 'homo_sapiens', 'nxnl2'
    # TODO, take care of aliases for species names, always use the
    # binomial names when running the code

    # 2-
    print("Searching ID for gene with name %s in species %s..." % (symbol,
                                                                   species))
    a = get_geneids_from_symbol(species, symbol)
    print("Found the following list of ids: %s" % (json.dumps(a)))
    if len(a) == 0:
        raise KeyError("No gene found, exiting")
    curgene = a[0]
    gene_name = "%s_%s" % (symbol, curgene)
    cdirectory = gene_name
    TSLsubdir = cdirectory + "/TSL/"
    TEXsubdir = cdirectory + "/TablesExons/"
    Seqsubdir = cdirectory + "/Sequences/"
    print("Using gene id %s from now on." % (curgene))
    print("Results will be saved in directory %s" % (cdirectory))
    if not os.path.exists(cdirectory):
        os.makedirs(TSLsubdir)
        os.makedirs(TEXsubdir)
        os.makedirs(Seqsubdir)
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
    nt = 0
    corthologs = Counter([x['target']['species'] for x in orthologs])
    for tl, c in corthologs.most_common():
        print("  %-23s: %4d" % (tl, c))
        # if nt > 5: break
        nt += 1
    ##
    orthologs_filtered = filter_ortho(
        orthokeep,
        orthologs,
        relation=["ortholog_one2one", "ortholog_one2many"])
    print("Filtering on %d species, %d matches" % (len(orthokeep),
                                                   len(orthologs_filtered)))
    ##
    print("Getting all the transcripts for preparing a TSL file")
    TSL_cur, TSL_ortho = get_transcripts_orthologs(curgene, orthologs_filtered)
    TSL_out = "%s/%s" % (TSLsubdir, gene_name)
    write_TSL_file(TSL_out, [TSL_cur] + TSL_ortho)

    print("**** Query species : %s" % (species))
    print("Got a total of %d transcripts with biotypes" % (len(TSL_cur)))
    for b, c in Counter([x['biotype'] for x in TSL_cur]).most_common():
        print("  %-23s: %4d" % (b, c))
    print("**** Orthologues")
    for tr_o in TSL_ortho:
        print("%-22s: %4d transcripts " % (tr_o[0]['species'], len(tr_o)))

    print("Now getting the exons sequences")
    # TODO revert to multiple files if it is easier
    ffasta = "%s/%s.fasta" % (Seqsubdir, gene_name)
    fexonstable = "%s/%s_exonstable.tsv" % (TEXsubdir, gene_name)
    fastaout = open(ffasta, "w")
    exonstableout = open(fexonstable, "w")
    dex = get_listofexons(curgene)
    lexid = list({x['exon_id'] for x in dex})
    print("Getting the sequences files for %s" % (curgene))
    curspec_ens_dataset = species2ensembldataset(species)
    exfasta = get_exons_sequences(lexid)
    extable = get_biomart_exons_annot(curspec_ens_dataset, curgene)
    exonstableout.write(extable.text)
    exons_name = "%s:%s" % (species, symbol)
    for dseq in exfasta:
        dictseq2fasta(dseq, exons_name, fastaout)
    for o in orthologs_filtered:
        orthoid = o['target']['id']
        orthospecies = o['target']['species']
        # orthotaxon = o['target']['taxon_id']
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
