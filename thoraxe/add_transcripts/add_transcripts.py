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
#from thoraxe.transcript_info import read_exon_file
from thoraxe.version import __version__


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
            raise Exception(f'The {col} column is missing in the input table.')
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
                f'Error with {name}: Species name should be lowercase, {example}'
            )
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


def calculate_cDNA_start_end(df):
    # Sort the dataframe by TranscriptID and ExonRank
    df = df.sort_values(by=['TranscriptID', 'ExonRank']).copy()
    
    # Convert necessary columns to numeric type
    numeric_columns = ['ExonRank', 'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart', 'GenomicCodingEnd', 'Strand']
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Initialize cDNA coding start and end columns
    df['cDNA_CodingStart'] = np.nan
    df['cDNA_CodingEnd'] = np.nan

    # Process each transcript separately
    for transcript_id, transcript_group in df.groupby('TranscriptID'):
        transcript_group = transcript_group.sort_values(by='ExonRank')

        # Compute exon lengths
        exon_lengths = (transcript_group['ExonRegionEnd'] - transcript_group['ExonRegionStart'] + 1).values
        
        # Identify coding exons
        is_coding_exon = transcript_group['GenomicCodingStart'].notna().values
        
        # Compute coding exon lengths
        coding_lengths = np.zeros_like(exon_lengths)
        coding_lengths[is_coding_exon] = (
            transcript_group.loc[is_coding_exon, 'GenomicCodingEnd'] - 
            transcript_group.loc[is_coding_exon, 'GenomicCodingStart'] + 1
        ).values
        
        # Compute 3' non-coding region
        nc_3_p = np.zeros_like(exon_lengths)
        strand = transcript_group['Strand'].iloc[0]  # All exons share the same strand
        
        if strand == -1:
            nc_3_p[is_coding_exon] = (
                transcript_group.loc[is_coding_exon, 'ExonRegionStart'] - 
                transcript_group.loc[is_coding_exon, 'GenomicCodingStart']
            ).values
        else:
            nc_3_p[is_coding_exon] = (
                transcript_group.loc[is_coding_exon, 'ExonRegionEnd'] - 
                transcript_group.loc[is_coding_exon, 'GenomicCodingEnd']
            ).values
        
        # Compute cumulative exon length
        cum_exon_len = np.cumsum(exon_lengths)

        # Compute cDNA coding start and end positions
        cDNA_coding_start = np.full_like(exon_lengths, np.nan, dtype=float)
        cDNA_coding_end = np.full_like(exon_lengths, np.nan, dtype=float)

        cDNA_coding_start[is_coding_exon] = (cum_exon_len - coding_lengths - nc_3_p + 1)[is_coding_exon]
        cDNA_coding_end[is_coding_exon] = (cum_exon_len - nc_3_p)[is_coding_exon]

        # Assign values back to the dataframe
        df.loc[transcript_group.index, 'cDNA_CodingStart'] = cDNA_coding_start
        df.loc[transcript_group.index, 'cDNA_CodingEnd'] = cDNA_coding_end

    return df

def clean_dataframe(input_df):
    """
    Cleans and prepares input DataFrame:
    - Calculates cDNA start and end values
    - Removes NaN values
    - Converts numeric columns to Int64
    - Sorts by TranscriptID and ExonRank
    """
    # First, calculate cDNA start and end values
    input_df = calculate_cDNA_start_end(input_df)
    
    # Remove NaN values and convert numeric columns
    #input_df = input_df.dropna().copy()
    numeric_columns = [
        'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart', 
        'GenomicCodingEnd', 'StartPhase', 'EndPhase', 'cDNA_CodingStart', 'cDNA_CodingEnd'
    ]
    
    for col in numeric_columns:
        if col in input_df.columns:
            input_df[col] = pd.to_numeric(input_df[col], errors='coerce').astype('Int64')
    
    # Sort data
    input_df = input_df.sort_values(by=['TranscriptID', 'ExonRank']).reset_index(drop=True)
    
    return input_df

def merge_with_exontable(input_df, exontable):
    """
    Merges cleaned input_df with exontable.
    """
    input_df = input_df.sort_values(by=['TranscriptID', 'ExonRank']).reset_index(drop=True)

    for row in input_df.itertuples(index=False):  # Dodato index=False
        exontable = pd.concat([
            exontable,
            pd.DataFrame({
                'GeneID': [row.GeneID],
                'TranscriptID': [row.TranscriptID],
                'ProteinID': [f"{row.TranscriptID}_PROTEIN"],
                'Strand': [row.Strand],
                'ExonID': [row.ExonID],
                'ExonRegionStart': [row.ExonRegionStart],
                'ExonRegionEnd': [row.ExonRegionEnd],
                'ExonRank': [row.ExonRank],
                'cDNA_CodingStart': [getattr(row, 'cDNA_CodingStart', None)],  # Korišćenje getattr()
                'cDNA_CodingEnd': [getattr(row, 'cDNA_CodingEnd', None)],
                'GenomicCodingStart': [getattr(row, 'GenomicCodingStart', None)],
                'GenomicCodingEnd': [getattr(row, 'GenomicCodingEnd', None)],
                'StartPhase': [getattr(row, 'StartPhase', None)],
                'EndPhase': [getattr(row, 'EndPhase', None)],
            })
        ], ignore_index=True)

    return exontable


def read_exon_file(exon_table_file):
    """Read the exon_table_file and return a cleaned pandas DataFrame."""
    
    int_cols_with_nas = [
        'cDNA_CodingStart', 'cDNA_CodingEnd', 'GenomicCodingStart',
        'GenomicCodingEnd'
    ]
    
    # Read CSV with dtype=str to avoid unwanted float conversion
    exon_data = pd.read_csv(exon_table_file, sep='\t', dtype='string')

    # Convert numeric columns to proper integer format, handling NaNs
    for col in int_cols_with_nas:
        if col in exon_data.columns:
            exon_data[col] = pd.to_numeric(exon_data[col], errors='coerce').astype('Int64')

    # Sort exon by rank in transcript
    exon_data = exon_data.sort_values(by=['GeneID', 'TranscriptID', 'ExonRank']).reset_index(drop=True)
    
    return exon_data

def read_transcript_file(transcript_file):
    """
    Read the transcript CSV file and return a pandas DataFrame without modifications.
    """
    # Read data from CSV
    return pd.read_csv(transcript_file)


def add_to_tsl(input_df, tsl_df):
    """
    It adds the new transcripts to the TSL table downloaded from Ensembl.
    """
    subset = input_df[['Species', 'TranscriptID']].drop_duplicates()

    # Create a new DataFrame to store the new transcript entries
    new_entries = pd.DataFrame({
        'Species': subset['Species'],
        'Name': subset['TranscriptID'],
        'TranscriptID': subset['TranscriptID'],
        'Source': 'user',
        'ExperimentSource': 'user',
        'Biotype': 'protein_coding',
        'Flags': np.nan,
        'Version': datetime.datetime.now().isoformat(timespec='hours')
    })

    # Use pd.concat() correctly
    tsl_df = pd.concat([tsl_df, new_entries], ignore_index=True)

    return tsl_df


def add_sequences(input_df, seqrecords):
    """
    It modifies seqrecords by appending the new sequences.
    """
    for row in input_df.itertuples():
        name = f'{row.Species}:{row.GeneID}'
        seqrecords.append(
            SeqRecord(
                Seq(row.NucleotideSequence, SingleLetterAlphabet()),
                id=name,
                name=name,
                description=
                f'{row.Species}:{row.GeneID} {row.ExonID} na:na:na:{row.ExonRegionStart}:{row.ExonRegionEnd}:{row.Strand}' # pylint: disable=line-too-long
            ))

def main():
    """Main script function to add user transcript data."""
    args = parse_command_line().parse_args()
    paths = get_output_paths(args)

    print("Reading user input...")
    input_df = read_transcript_file(paths.input)  

    exontable = read_exon_file(paths.exontable)
    tsl_df = pd.read_csv(paths.tsl)

    print("Checking input...")
    check_input(input_df, exontable)

    input_df = clean_dataframe(input_df)

    print("Adding new transcripts...")
    new_exontable = merge_with_exontable(input_df, exontable)

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