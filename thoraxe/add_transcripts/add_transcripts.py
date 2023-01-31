"""
Functions to add user-defined transcript data from a CSV file to the
previously downloaded Ensembl data.
"""

import argparse
import collections
import datetime
import logging
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

from tabulate import tabulate
from thoraxe.transcript_info import read_exon_file
from ..version import __version__


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse add_transcripts' command line arguments and
    returns the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="add_transcripts",
        description="""
        add_transcripts downloads add user-defined transcripts to
        the data previously download from Ensembl.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input',
                        help='Input CSV containing the transcript data.',
                        type=str)
    parser.add_argument(
        'ensembl',
        help='Path to the previously download Ensembl data for a gene, '
        'i.e. the path to the Ensembl directory created by '
        'transcript_query',
        type=str)
    parser.add_argument('-v',
                        '--verbose',
                        help='Print detailed progress.',
                        action='store_true')
    parser.add_argument('--version',
                        action='version',
                        version=f'ThorAxe version {__version__}')

    return parser


COLUMN_NAMES = {
    'Species',
    'GeneID',
    'TranscriptID',
    'Strand',
    'ExonID',
    'ExonRank',
    'ExonRegionStart',
    'ExonRegionEnd',
    'GenomicCodingStart',
    'GenomicCodingEnd',
    'StartPhase',
    'EndPhase',
    'NucleotideSequence',
}

Paths = collections.namedtuple('Paths',
                               ['input', 'tsl', 'exontable', 'sequences'])


def get_output_paths(args):
    """
    Return a `Paths` `namedtuple` containing the input and output paths for
    `add_transcripts`.
    """
    return Paths(args.input, os.path.join(args.ensembl, 'tsl.csv'),
                 os.path.join(args.ensembl, 'exonstable.tsv'),
                 os.path.join(args.ensembl, 'sequences.fasta'))


def read_fasta(sequences_path):
    """
    Returns a list of `SeqRecord`s for the fasta file at `sequences_path`.
    """
    return list(SeqIO.parse(sequences_path, "fasta"))


def _write_fasta(sequences_path, seqrecords):
    SeqIO.write(seqrecords, sequences_path, "fasta")


def _get_exon_sequences(seqrecords):
    return {
        seqrecord.description.split(' ')[1]: str(seqrecord.seq)
        for seqrecord in seqrecords
    }


def _in(value, column):
    return any(value == column)


def _check_col_names(input_df):
    for col in COLUMN_NAMES:
        if col not in input_df.columns:
            raise Exception(
                f'The {col} column is missing in the input table.')
    for col in input_df.columns:
        if col not in COLUMN_NAMES:
            logging.warning('The %s column is not going to be used', col)


def _check_transcript(input_df, exontable):
    for transcript in input_df['TranscriptID'].unique():
        if _in(transcript, exontable['TranscriptID']):
            raise Exception(
                f'TranscriptID {transcript} is already in the Ensembl data!')


def _format_df(data):
    return tabulate(data,
                    headers='keys',
                    tablefmt='psql',
                    floatfmt=".0f",
                    showindex=False)


def _check_exon(input_df, exontable):
    exon_columns = [
        'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart',
        'GenomicCodingEnd', 'StartPhase', 'EndPhase'
    ]
    for exon in input_df['ExonID'].unique():
        if _in(exon, exontable['ExonID']):
            logging.warning(
                "%s is on the Ensembl data; the new sequence won't be used.",
                exon)
            previous_exon = exontable[exontable.ExonID ==
                                      exon][exon_columns].drop_duplicates()
            new_exon = input_df[input_df.ExonID ==
                                exon][exon_columns].drop_duplicates()
            if not (previous_exon.values == new_exon.values).all():
                raise Exception(f'''
Exon {exon} is already in the Ensembl data with different values!
                    
Exon data at Ensembl: 
{_format_df(previous_exon)}

New exon data:
{_format_df(new_exon)}''')


def _check_species_name(input_df):
    for row in input_df.itertuples():
        example = 'e.g. homo_sapiens.'
        name = row.Species
        if not name.islower():
            raise Exception(
                f'Error with {name}: Species name should be lowercase, {example}')
        if len(name.split('_')) != 2:
            raise Exception(
                f'Error with {name}:' +
                'Species name should be binomial, and the terms should be '
                f'separated by underscore, {example}')


def _check_pair_order(input_df):
    for row in input_df.itertuples():
        if row.ExonRegionStart > row.ExonRegionEnd:
            raise Exception(
                'ExonRegionStart should be lower than ExonRegionEnd.')
        if row.GenomicCodingStart > row.GenomicCodingEnd:
            raise Exception(
                'GenomicCodingStart should be lower than GenomicCodingEnd.')
        if len(row.NucleotideSequence) != (row.ExonRegionEnd -
                                           row.ExonRegionStart + 1):
            raise Exception(
                "The sequence length estimated using ExonRegionStart and "
                "ExonRegionEnd disagrees with the actual NucleotideSequence "
                "length.")


def check_input(input_df, exontable):
    """
    It raises an error if there is a problem with the user-defined transcripts.
    """
    _check_col_names(input_df)
    _check_species_name(input_df)
    _check_transcript(input_df, exontable)
    _check_pair_order(input_df)
    _check_exon(input_df, exontable)


def add_to_exontable(input_df, exontable):
    """
    It adds the new exons and transcripts to the Ensembl's exontable.
    """
    transcript_id = ""
    cdna_start = 0
    cdna_end = 0
    for row in input_df.itertuples():
        if transcript_id != row.TranscriptID or cdna_end == 0:
            transcript_id = row.TranscriptID
            coding_len = (row.GenomicCodingEnd - row.GenomicCodingStart)
            exon_len = (row.ExonRegionEnd - row.ExonRegionStart)
            cdna_start = exon_len - coding_len + 1
            cdna_end = cdna_start + coding_len
        else:
            cdna_start = cdna_end + 1
            cdna_end += (row.GenomicCodingEnd - row.GenomicCodingStart + 1)
        exontable = exontable.append(
            {
                'GeneID': row.GeneID,
                'TranscriptID': row.TranscriptID,
                'ProteinID': f"{row.TranscriptID}_PROTEIN",
                'Strand': row.Strand,
                'ExonID': row.ExonID,
                'ExonRegionStart': row.ExonRegionStart,
                'ExonRegionEnd': row.ExonRegionEnd,
                'ExonRank': row.ExonRank,
                'cDNA_CodingStart': cdna_start,
                'cDNA_CodingEnd': cdna_end,
                'GenomicCodingStart': row.GenomicCodingStart,
                'GenomicCodingEnd': row.GenomicCodingEnd,
                'StartPhase': row.StartPhase,
                'EndPhase': row.EndPhase,
            },
            ignore_index=True)
    return exontable


def add_to_tsl(input_df, tsl_df):
    """
    It adds the new transcripts to the TSL table downloaded from Ensembl.
    """
    subset = input_df[['Species', 'TranscriptID']].drop_duplicates()
    for row in subset.itertuples():
        tsl_df = tsl_df.append(
            {
                'Species': row.Species,
                'Name': row.TranscriptID,
                'TranscriptID': row.TranscriptID,
                'Source': 'user',
                'ExperimentSource': 'user',
                'Biotype': 'protein_coding',
                'Flags': np.nan,
                'Version': str(datetime.datetime.now().isoformat(timespec='hours')),
            },
            ignore_index=True)
    return tsl_df


def add_sequences(input_df, seqrecords):
    """
    It modifies seqrecords by appending the new sequences.
    """
    for row in input_df.itertuples():
        name = f'{row.Species}:{row.GeneID}'
        seqrecords.append(
            SeqRecord(Seq(row.NucleotideSequence, SingleLetterAlphabet()),
                      id=name,
                      name=name,
                      description=f'{row.Species}:{row.GeneID} {row.ExonID} na:na:na:{row.ExonRegionStart}:{row.ExonRegionEnd}:{row.Strand}'))


def main():
    """Main script function to add user transcript data."""
    args = parse_command_line().parse_args()
    paths = get_output_paths(args)
    print("Reading user input...")
    input_df = pd.read_csv(paths.input)
    exontable = read_exon_file(paths.exontable)
    tsl_df = pd.read_csv(paths.tsl)
    print("Checking input...")
    check_input(input_df, exontable)
    print("Adding new transcripts...")
    new_exontable = add_to_exontable(input_df, exontable)
    new_tsl = add_to_tsl(input_df, tsl_df)
    seqrecords = read_fasta(paths.sequences)
    add_sequences(input_df, seqrecords)
    print("Saving...")
    _write_fasta(paths.sequences, seqrecords)
    new_tsl.to_csv(paths.tsl, index=False)
    new_exontable.to_csv(paths.exontable, sep='\t', index=False)
    print("Finished!")


if __name__ == '__main__':
    main()
