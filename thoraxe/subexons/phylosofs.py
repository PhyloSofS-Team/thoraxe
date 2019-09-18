"""
phylosofs: Functions to generate PhyloSofS input.
"""

import os
import string
import warnings

import pandas as pd

from Bio import Phylo
from Bio import BiopythonWarning
from thoraxe.utils import folders
from thoraxe.utils import split

CHARS = [
    char for char in string.printable if char not in {
        ' ', '\t', '\n', '\r', '\x0b', '\x0c', '\\', '*', '>', '"', "'", ',',
        '-', '_', '/', ';', '#', '$', '.', '&', '!', '@'
    }
]

CHARS.extend([
    '¼', '½', '¾', '¿', 'À', 'Á', 'Â', 'Ã', 'Ä', 'Å', 'Æ', 'Ç', 'È', 'É', 'Ê',
    'Ë', 'Ì', 'Í', 'Î', 'Ï', 'Ð', 'Ñ', 'Ò', 'Ó', 'Ô', 'Õ', 'Ö', '×', 'Ø', 'Ù',
    'Ú', 'Û', 'Ü', 'Ý', 'Þ', 'ß', 'à', 'á', 'â', 'ã', 'ä', 'å', 'æ', 'ç', 'è',
    'é', 'ê', 'ë', 'ì', 'í', 'î', 'ï', 'ð', 'ñ', 'ò', 'ó', 'ô', 'õ', 'ö', '÷',
    'ø', 'ù', 'ú', 'û', 'ü', 'ý', 'þ', 'ÿ', 'Ā', 'ā', 'Ă', 'ă', 'Ą', 'ą', 'Ć',
    'ć', 'Ĉ', 'ĉ', 'Ċ', 'ċ', 'Č', 'č', 'Ď', 'ď', 'Đ', 'đ', 'Ē', 'ē', 'Ĕ', 'ĕ',
    'Ė', 'ė', 'Ę', 'ę', 'Ě', 'ě', 'Ĝ', 'ĝ', 'Ğ', 'ğ', 'Ġ', 'ġ', 'Ģ', 'ģ', 'Ĥ',
    'ĥ', 'Ħ', 'ħ', 'Ĩ', 'ĩ', 'Ī', 'ī', 'Ĭ', 'ĭ', 'Į', 'į', 'İ', 'ı', 'Ĳ', 'ĳ',
    'Ĵ', 'ĵ', 'Ķ', 'ķ', 'ĸ', 'Ĺ', 'ĺ', 'Ļ', 'ļ', 'Ľ', 'ľ', 'Ŀ', 'ŀ', 'Ł', 'ł',
    'Ń', 'ń', 'Ņ', 'ņ', 'Ň', 'ň', 'ŉ', 'Ŋ', 'ŋ', 'Ō', 'ō', 'Ŏ', 'ŏ', 'Ő', 'ő',
    'Œ', 'œ', 'Ŕ', 'ŕ', 'Ŗ', 'ŗ', 'Ř', 'ř', 'Ś', 'ś', 'Ŝ', 'ŝ', 'Ş', 'ş', 'Š',
    'š', 'Ţ', 'ţ', 'Ť', 'ť', 'Ŧ', 'ŧ', 'Ũ', 'ũ', 'Ū', 'ū', 'Ŭ', 'ŭ', 'Ů', 'ů',
    'Ű', 'ű', 'Ų', 'ų', 'Ŵ', 'ŵ', 'Ŷ', 'ŷ', 'Ÿ', 'Ź', 'ź', 'Ż', 'ż', 'Ž', 'ž',
    'ſ', 'Ơ', 'ơ', 'Ư', 'ư', 'Ǆ', 'ǅ', 'ǆ', 'Ǉ', 'ǈ', 'ǉ', 'Ǌ', 'ǋ', 'ǌ', 'Ǎ',
    'ǎ', 'Ǐ', 'ǐ', 'Ǒ', 'ǒ', 'Ǔ', 'ǔ', 'Ǖ', 'ǖ', 'Ǘ', 'ǘ', 'Ǚ', 'ǚ', 'Ǜ', 'ǜ',
    'ǰ', 'Ǳ', 'ǲ', 'ǳ', 'Ǵ', 'ǵ', 'Ƕ', 'Ƿ', 'Ǹ', 'ǹ', 'Ǻ', 'ǻ', 'Ǽ', 'ǽ', 'Ǿ',
    'ǿ', 'Ȁ', 'ȁ', 'Ȃ', 'ȃ', 'Ȅ', 'ȅ', 'Ȇ', 'ȇ', 'Ȉ', 'ȉ', 'Ȋ', 'ȋ', 'Ȍ', 'ȍ',
    'Ȏ', 'ȏ', 'Ȑ', 'ȑ', 'Ȓ', 'ȓ', 'Ȕ', 'ȕ', 'Ȗ', 'ȗ', 'Ș', 'ș', 'Ț', 'ț', 'Γ',
    'Δ', 'Θ', 'Λ', 'Ξ', 'Π', 'Σ', 'Φ', 'Ψ', 'Ω', 'Ϊ', 'Ϋ', 'ά', 'έ', 'ή', 'ί',
    'ΰ', 'α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ', 'λ', 'μ', 'ν', 'ξ',
    'ο', 'π', 'ρ', 'ς', 'σ', 'τ', 'υ', 'φ', 'χ', 'ψ', 'ω', 'ϊ', 'ϋ', 'ό', 'ύ',
    'ώ', 'ϑ', 'ϒ', 'ϓ', 'ϖ', 'Ϣ', 'ϣ', 'Ϥ', 'ϥ', 'Ϧ', 'ϧ', 'Ϩ'
])

CHARS.extend([
    '㐀', '㐁', '㐂', '㐃', '㐄', '㐅', '㐆', '㐇', '㐈', '㐉', '㐊', '㐧', '㐪', '㐫', '㐬',
    '㐭', '㐱', '㐴', '㐵', '㐶', '㞥', '㞦', '㞧', '㞨', '㞩', '㞪', '㞫', '㞬', '㞭', '䇂'
])

CHARS = sorted(set(CHARS))

MAX_EXONS = len(CHARS)


def get_exon2char(exons):
    """
    Return a dictionary from orthologous exon group to a single character.

    >>> result = get_exon2char(['1_0/1_1', '2_0'])
    >>> sorted(result)
    ['1_0', '1_1', '2_0']
    >>> [result[exon] for exon in sorted(result)]
    ['%', '(', ')']
    """
    homologous_exons = {}
    i = 0
    for exon in exons:
        for subexon in str(exon).split('/'):
            if subexon not in homologous_exons:
                if i < MAX_EXONS:
                    homologous_exons[subexon] = CHARS[i]
                    i += 1
                else:
                    raise Exception(
                        'PhyloSofS can parse more than {} homologous exons.'.
                        format(MAX_EXONS))
    return homologous_exons


def get_transcript2phylosofs(data, exon2char):
    """
    Return a dictionary from transcript id to its phylosofs representation.
    """
    transcript2phylosofs = {}
    for _, gene in data.groupby('GeneID'):
        for transcript_id, transcript in gene.groupby('TranscriptID'):
            phylosofs = []
            for subexon in transcript.sort_values(
                    'SubexonRank')['HomologousExons']:
                for ortholog in subexon.split('/'):
                    phylosofs.append(exon2char[ortholog])
            transcript2phylosofs[transcript_id] = ''.join(phylosofs)
    return transcript2phylosofs


def _get_terminal_names(input_tree):
    """Return the set of terminal names in the Newick tree at input_tree."""
    tree = Phylo.read(input_tree, "newick")
    return {clade.name for clade in tree.get_terminals()}


def _get_protein2gene(exontable_file):
    """Return a pandas datafarme to map protein (translation) id to gene id."""
    data = pd.read_csv(exontable_file, sep='\t')
    with warnings.catch_warnings():
        # Bio/Seq.py : class Seq : __hash__ : warnings.warn
        warnings.simplefilter('ignore', BiopythonWarning)
        data = data.loc[:, ['GeneID', 'ProteinID']].dropna().drop_duplicates(
        ).set_index('ProteinID')

    data['GeneID'] = data['GeneID'].astype('category')
    return data


def _get_terminals_to_delete(tree, used_proteins):
    """Return the list of terminal nodes to delete."""
    return [
        clade.name for clade in tree.get_terminals()
        if clade.name not in used_proteins and clade.name is not None
    ]


def prune_tree(input_tree, output_tree, exontable_file, used_genes):
    """
    Delete unused proteins from the Ensembl gene tree.

    The function reads the Newick tree from input_tree and stores the pruned
    tree in output_tree and returning the output_tree path. If an
    error occurs, the function does not save the tree and raises an Exception.

    In the tree, each protein is identified by their ensembl translation id.
    For each gene, Ensembl uses the longest translated transcript sequence.

    The exontable_file downloaded from Ensembl with transcript_query is needed
    for getting the mapping between translation and gene ids of the used_genes.
    """
    tree = Phylo.read(input_tree, "newick")
    data = _get_protein2gene(exontable_file)

    # Delete unused proteins
    used_proteins = set(data.index[data['GeneID'].isin(used_genes)])
    for clade in _get_terminals_to_delete(tree, used_proteins):
        try:
            tree.prune(clade)
        except ValueError as err:
            raise Exception(
                'Error with protein {} while prunning the tree {} : {}'.format(
                    clade, tree, err))

    # Use gene id instead of protein/peptide/translation id in the tree
    for clade in tree.get_terminals():
        clade.name = data.loc[clade.name]

    Phylo.write(tree, output_tree, 'newick')

    return output_tree


def _fill_sequence_and_annotation(df_group, exon2char):
    """Create a list of sequences and homologous exons (annotation)."""
    exon_annot = []
    seqs = []
    for row in df_group.itertuples():
        seq = str(row.SubexonProteinSequence).replace('*', '')
        exons = split.split_exons(row.HomologousExons)
        exon_lengths = split.split_lengths(row.HomologousExonLengths, seq)
        for exon, length in zip(exons, exon_lengths):
            for _ in range(length):
                exon_annot.append(exon2char[exon])
        seqs.append(seq)
    return "".join(seqs), "".join(exon_annot)


def _transcript_pir(exon_data, output_file, exon2char, transcript2phylosofs):
    """
    Create a PIR file with transcript sequences.

    Annotation line is used to store the one character identifier of the
    homologous exon for each residue.
    """
    with open(output_file, 'w', encoding='utf-8') as file:
        with warnings.catch_warnings():
            # Bio/Seq.py : class Seq : __hash__ : warnings.warn
            warnings.simplefilter('ignore', BiopythonWarning)
            groups = exon_data.loc[:, [
                'GeneID', 'TranscriptID', 'SubexonProteinSequence',
                'SubexonRank', 'HomologousExons', 'HomologousExonLengths',
                'HomologousExonSequences'
            ]].drop_duplicates().sort_values(
                by=['GeneID', 'TranscriptID', 'SubexonRank']).groupby(
                    ['GeneID', 'TranscriptID'])
        for (gene, transcript), group in groups:
            seq, annot = _fill_sequence_and_annotation(group, exon2char)
            file.write(">P1;{} {} {}\n".format(
                gene, transcript, transcript2phylosofs[transcript]))
            file.write(annot + "\n")
            file.write(seq + "*\n")
    return output_file


def _save_transcripts(exon_data, output_file, transcript2phylosofs):
    """Save transcripts in PhyloSofS' format."""
    with open(output_file, 'w', encoding='utf-8') as file:
        gene_transcripts = exon_data.loc[:, ['GeneID', 'TranscriptID'
                                             ]].sort_values(
                                                 by=['GeneID', 'TranscriptID'
                                                     ]).drop_duplicates()
        for gene, group in gene_transcripts.groupby('GeneID'):
            file.write("{} ".format(gene))
            transcripts = group['TranscriptID']
            n_transcripts = len(transcripts)
            for i in range(n_transcripts):
                file.write(transcript2phylosofs[transcripts.iloc[i]])
                if i != n_transcripts - 1:
                    file.write(" ")
                else:
                    file.write("\n")
    return output_file


def phylosofs_inputs(exon_data, ensembl_folder, output_folder):
    """
    Create the input files needed for running PhyloSofS in a phylosofs folder.
    """
    output_path = folders.create_subfolder(output_folder, 'phylosofs')

    exon2char = get_exon2char(exon_data['HomologousExons'].unique())
    transcript2phylosofs = get_transcript2phylosofs(exon_data, exon2char)

    prune_tree(os.path.join(ensembl_folder, 'tree.nh'),
               os.path.join(output_path, 'tree.nh'),
               os.path.join(ensembl_folder, 'exonstable.tsv'),
               list(exon_data['GeneID'].unique()))

    with open(os.path.join(output_path, 'homologous_exons.tsv'),
              'w',
              encoding='utf-8') as file:
        for exon, char in exon2char.items():
            file.write("{}\t{}\n".format(exon, char))

    _save_transcripts(exon_data, os.path.join(output_path, 'transcripts.txt'),
                      transcript2phylosofs)

    _transcript_pir(exon_data, os.path.join(output_path, 'transcripts.pir'),
                    exon2char, transcript2phylosofs)
