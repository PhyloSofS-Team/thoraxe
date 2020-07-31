"""Functions that use pandas to read transcript information."""

import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import pandas as pd

from thoraxe.transcript_info import clusters
from thoraxe.transcript_info import phases


def _check_column_presence(data_frame, column, message=""):
    """Raise an error if column is not in data_frame."""
    if column not in data_frame.columns:
        raise ValueError("Input DataFrame hasn't column: %s. %s" %
                         (column, message))


def _check_column_absence(data_frame, column, message=""):
    """Give a warning if the column is already present in the data_frame."""
    if column in data_frame.columns:
        logging.warn("Input DataFrame already has column: %s. %s", column,
                     message)


def _get_flag(flag):
    """
    Get TSL flag as float.

    It takes a Transcript Support Level (TSL) `flag`/value and returns it as
    a float if the input is a str.

    >>> _get_flag("5")
    5.0
    >>> _get_flag(5.0)
    5.0
    >>> _get_flag("1 (assigned to previous version 6)")
    1.0
    >>> import numpy as np
    >>> np.isnan(_get_flag("NA"))
    True
    """
    if isinstance(flag, str):
        flag_code = flag.split(' ')[0]
        if flag_code.upper() in ['NA', 'NAN']:
            return np.nan
        return float(flag_code)

    return flag


def _in_species_list(species, species_list=None):
    """
    Return True if species_list is None or if species is in species_list.

    >>> _in_species_list('homo_sapiens')
    True
    >>> _in_species_list('homo_sapiens', ['homo_sapiens', 'mus_musculus'])
    True
    >>> _in_species_list('pan_paniscus', ['homo_sapiens', 'mus_musculus'])
    False
    """
    if species_list is None:
        return True
    return species in species_list


def read_tsl_file(tsl_file, max_tsl_level, remove_na=False, species_list=None):
    """
    Read a csv file with the Transcript Support Level (TSL) data.

    `max_tsl_level` determines the maximum Transcript Support Level (TSL)
    level to keep. The options are:

    - 1.0 : All splice junctions of the transcript are supported by at
      least one non-suspect mRNA.
    - 2.0 : The best supporting mRNA is flagged as suspect or the support is
      from multiple ESTs (Expressed Sequence Tags).
    - 3.0 : The only support is from a single EST.
    - 4.0 : The best supporting EST is flagged as suspect.
    - 5.0 : No single transcript supports the model structure.

    If `remove_na` is True, columns with missing TSL are not included
    (default: False).

    It returns a pandas DataFrame with MultiIndex: ('Species', 'TranscriptID').
    """
    assert max_tsl_level in {1.0, 2.0, 3.0, 4.0, 5.0}

    # TO DO : Modify EnsemblRESTTranscriptQueries.py to make stable the input
    #         of this function:
    # TO DO :     1 - To ensure column names,
    # TO DO :     2 - To ensure that Biotype always use 'protein_coding' to
    #                 simplify the selection.

    tsl_data = pd.read_csv(tsl_file)

    # Select only protein coding transcrits :
    # ---------------------------------------
    # QUESTION : Can Biotype be 'Protein coding' or any other thing rather
    # than 'protein_coding'?
    # ANSWER : Until now we only saw 'Protein coding' or 'protein_coding'
    # in the ENSEMBL annotation but there will be surely other variations
    # on those words.

    tsl_data = tsl_data[[
        _in_species_list(row.Species, species_list)
        and row.Biotype in ['protein_coding', 'Protein coding']
        for row in tsl_data.itertuples()
    ]]

    # explicit convert to float, sometimes Flags are read as str
    tsl_data["TSL"] = tsl_data.Flags.map(_get_flag)
    tsl_data.drop(columns=['Flags'], inplace=True)
    selected_flags = tsl_data.TSL <= max_tsl_level

    if remove_na:
        to_delete = set(tsl_data.TranscriptID[~selected_flags])
        tsl_data = tsl_data[selected_flags]
    else:
        missing_tsl = tsl_data.TSL.isna()
        to_delete = set(tsl_data.TranscriptID[~(selected_flags | missing_tsl)])
        tsl_data = tsl_data[selected_flags | missing_tsl]

    if to_delete:
        clusters.inform_about_deletions(
            to_delete,
            "There are transcripts with TSL greater than {}{}:".format(
                max_tsl_level, " or no TSL information" if remove_na else ""))

    return tsl_data


def read_exon_file(exon_table_file):
    """Read the exon_table_file and returns a pandas' DataFrame."""
    # Read data :
    # -----------
    int_cols_with_nas = [
        'cDNA_CodingStart', 'cDNA_CodingEnd', 'GenomicCodingStart',
        'GenomicCodingEnd'
    ]
    exon_data = pd.read_csv(exon_table_file,
                            sep='\t',
                            dtype={col: np.str
                                   for col in int_cols_with_nas})
    # That columns are sometimes interpreted as
    # floats because pandas doesn't allow NA values
    # for int columns. We read these columns as np.str
    # (object) to later delete rows with NAs and then
    # convert values to np.int.
    exon_data.dropna(inplace=True)

    for col in int_cols_with_nas:
        exon_data[col] = pd.to_numeric(exon_data[col])

    # Sort exon by rank in transcript :
    # ---------------------------------
    exon_data.sort_values(by=['GeneID', 'TranscriptID', 'ExonRank'],
                          inplace=True)

    return exon_data


def _complete_transcripts(data):
    """
    Delete transcripts that have exons without sequence.

    There are exons in the TSL table that doesn't have sequences in the
    downloaded fasta file, e.g.: ENSBBBE00000162100 for MAPK10
    """
    mask = data.ExonSequence.isnull()
    if mask.sum() > 0:
        exons_to_delete = set(data.ExonID[mask])
        to_delete = set(data.TranscriptID[mask])
        index_to_delete = data.index[data.TranscriptID.isin(to_delete)]
        data.drop(index_to_delete, inplace=True)
        clusters.inform_about_deletions(
            to_delete,
            "Exons without sequences: {}... Deleting incomplete transcripts:".
            format(exons_to_delete))


def add_exon_sequences(data_frame, sequence_file):
    """
    Add exon sequences to the DataFrame.

    This function adds an 'ExonSequence' column at the end of 'data_frame'.
    For each row, this column has the BioPython's SeqRecord for the exon
    sequence of the row. The sequence description in the sequence_file
    (in fasta format) should have the 'ExonID' as the second element
    if the description is split by ' '.ipt_in
    The data_frame should have an 'ExonID' column and should not have
    an 'ExonSequence' column.
    'ExonID' is going to be used to match the SeqRecord to the row.
    """
    _check_column_presence(data_frame, 'ExonID')
    _check_column_absence(data_frame, 'ExonSequence',
                          'Values are going to change.')

    data_frame['ExonSequence'] = None
    for seqrecord in SeqIO.parse(sequence_file, "fasta"):
        exon_id = seqrecord.description.split(' ')[1]
        # NOTE : This takes 180 ms where one explicit loop over the rows
        # takes 2.75 s :
        selected_rows = data_frame.index[data_frame['ExonID'] ==
                                         exon_id].tolist()
        # NOTE : Assign one value / row at each time is needed to avoid
        # shape errors :
        for row in selected_rows:
            # NOTE : at is needed to assign the sequence in-place
            data_frame.at[row, 'ExonSequence'] = seqrecord.seq

    # Delete transcripts without sequences
    _complete_transcripts(data_frame)


def _get_exon_cds(exon_cdna_seq, cdna_coding_start, cdna_coding_end,
                  first_exon, last_exon):
    """
    Return the CDS of the exon.

    This function takes the cDNA sequence of the exon and the
    coding start and end coordinates to cut the sequence and
    return the coding sequence (CDS) of the exon.
    """
    # ExonSequences in the table are cDNA sequences.
    # cDNA sequence can have UTR regions but they have not introns.
    # CDS hasn't UTRs nor introns :ipt_in
    #
    # gene : [..first.exon..]---intron---[...last.exon...]
    # cDNA : [.UTR..|ATG....]            [.....|...UTR...]
    # CDS  :        |ATG....]            [.....|
    #
    # For sequences in the strand -1 the table has the reverse complement.
    #
    # cDNA_CodingStart and end are relative to the cDNA of the first exon.
    # However, the first exon is not always coding,
    # i.e. it could not be in the table.
    #
    #                       cds_len
    #                      ╔═══════╗
    #
    # cDNA : [.....UTR.....|ATG....]
    #        ^             ^
    #        cdna_start    cdna_coding_start
    #
    #        ╚═════════════════════╝
    #                cdna_len
    #

    cdna_len = len(exon_cdna_seq)
    cds_len = cdna_coding_end - cdna_coding_start + 1

    if cds_len == cdna_len:  # if there isn't UTR region
        return exon_cdna_seq
    if first_exon:
        return exon_cdna_seq[-cds_len:]
    if last_exon:
        return exon_cdna_seq[:cds_len]

    raise ValueError('exon has {} bases instead of {}'.format(
        cdna_len, cds_len))


def _is_first_or_last_exon(row_list, row_index):
    """Determine if it is the first or the last coding exon."""
    start_exon = False
    end_exon = False

    n_rows = len(row_list)
    if row_index == 0:
        start_exon = True
    if row_index == (n_rows - 1):
        end_exon = True

    row = row_list[row_index]
    if row_index > 0:
        identical_to_previous = row_list[
            row_index - 1]['TranscriptID'] == row['TranscriptID']
        if not identical_to_previous:  # or (row['StartPhase'] == -1):
            start_exon = True
    if row_index < (n_rows - 1):
        identical_to_next = row_list[row_index +
                                     1]['TranscriptID'] == row['TranscriptID']
        if not identical_to_next:  # or (row['EndPhase'] == -1):
            end_exon = True

    return start_exon, end_exon


def _manage_start_phase_negative_strand(cds_seq, start_phase):
    """
    Cut the CDS sequence depending on the exon start phase.

    If the phase is different from 0 and gene in the negative strand, the
    codon shared between exons is assigned to the exon at the left. That means
    that our CDS sequence should not include the partial codon at the
    beginning. This function returns the cds_seq without the partial codon when
    the start_phase is 1 or 2.

    >>> _manage_start_phase_negative_strand('XXXYYY', 0)
    'XXXYYY'
    >>> _manage_start_phase_negative_strand('XXYYY', 1)
    'YYY'
    >>> _manage_start_phase_negative_strand('XYYY', 2)
    'YYY'
    """
    if cds_seq is None:
        return None

    if start_phase == 1:
        return cds_seq[2:]
    if start_phase == 2:
        return cds_seq[1:]

    return cds_seq  # do not change CDS sequence if exon start phase is 0 or -1


def _manage_end_phase_negative_strand(row_list, row_index, cds_seq, end_exon):
    """
    Complete the last codon of the CDS sequence depending on the end phase.

    If the phase is different from 0 and the gene in the negative strand, the
    codon shared between exons is assigned to the exon at the left (this exon).
    That means that our CDS sequence should complete its partial codon using
    bases from the beginning of the next exon. This function returns the
    cds_seq with the complete codon at the end when the end_phase is 1 or 2.
    If the last exon of the transcript has an incomplete codon at the end, it is
    deleted on the returned sequence.

    NOTE: _check_phases_by_position should be used before this function
    to ensure that there are not errors in the input data.
    """
    if cds_seq is None:
        return None
    end_phase = row_list[row_index]['EndPhase']
    if end_phase in {1, 2}:
        if end_exon:
            cds_seq = cds_seq[:-end_phase]  # delete incomplete codon
        else:
            next_exon_sequence = row_list[row_index + 1]['ExonSequence']
            n_bases = phases.bases_to_complete_previous_codon(end_phase)
            cds_seq = cds_seq + next_exon_sequence[0:n_bases]

    return cds_seq


def _manage_start_phase_positive_strand(row_list, row_index, cds_seq,
                                        start_exon):
    """
    Complete the first codon of the CDS sequence depending on the start phase.

    If the phase is different from 0 and the gene in the positive strand, the
    codon shared between exons is assigned to the exon at the right (this exon).
    That means that our CDS sequence should complete its partial codon using
    bases from the end of the previous exon. This function returns the
    cds_seq with the complete codon at the beginning when the end_phase is 1 or
    2. If the first exon of the transcript has an incomplete codon at the
    beginning, it is deleted on the returned sequence.

    >>> # Example based in ENSPSME00000216679 and ENSPSME00000216683
    >>> seq_one = 'AATTGATATATAAGGAAGTTATGGACTTGGAAGAGAGAACCAAGAATGGAGTCATACGG'
    >>> seq_one = seq_one + 'GGGCAGCCTTCTCCTTTAGGTTG'
    >>> seq_two = 'TCAAAAGACACATTCTGTAGTGTTGTAA'
    >>> row_one = {'StartPhase': 1, 'EndPhase': 2, 'ExonSequence': seq_one}
    >>> row_two = {'StartPhase': 2, 'EndPhase': 0, 'ExonSequence': seq_two}
    >>> row_list = [row_one, row_two]
    >>> _manage_start_phase_positive_strand(row_list, 1, seq_two, False)
    'TGTCAAAAGACACATTCTGTAGTGTTGTAA'
    """
    if cds_seq is None:
        return None
    start_phase = row_list[row_index]['StartPhase']
    if start_phase in {1, 2}:
        if start_exon:
            n_bases = phases.bases_to_complete_previous_codon(start_phase)
            cds_seq = cds_seq[n_bases:]  # delete incomplete codon
        else:
            n_bases = phases.bases_to_complete_next_codon(start_phase)
            if n_bases is None:
                return None
            seq = row_list[row_index - 1]['ExonSequence']
            seq = seq[-n_bases:]  # pylint: disable=invalid-unary-operand-type
            cds_seq = seq + cds_seq

    return cds_seq


def _manage_end_phase_positive_strand(cds_seq, end_phase):
    """
    Cut the CDS sequence depending on the exon end phase.

    If the phase is different from 0 and the gene in the positive strand,
    the codon shared between exons is assigned to the exon at the right.
    That means that our CDS sequence should not include the partial codon
    at the end. This function returns the cds_seq without the partial codon
    when the end_phase is 1 or 2.

    >>> _manage_end_phase_positive_strand('XXXY', 1)
    'XXX'
    >>> _manage_end_phase_positive_strand('XXXYY', 2)
    'XXX'
    >>> _manage_end_phase_positive_strand('XXXYYY', 0)
    'XXXYYY'
    """
    if cds_seq is None:
        return None
    if end_phase in {1, 2}:
        n_bases = phases.bases_to_complete_next_codon(end_phase)
        if n_bases is None:
            return None
        return cds_seq[:-n_bases]  # pylint: disable=invalid-unary-operand-type

    return cds_seq  # do not change CDS sequence if exon start phase is 0 or -1


def _check_phases_by_position(row, end_exon, allow_incomplete_cds):
    """
    Raise an error if the end phase don't make sense for the exon position.

    In particular, an exon end phase -1 is only allowed for the last exon.
    If end_exon is True, the exon has the last position in the transcript.
    If allow_incomplete_cds is False, an error is raised when the last exon
    finishes with an incomplete codon (end_phase 1 or 2).
    The pandas' DataFrame row is printed in the error message.
    """
    end_phase = row['EndPhase']
    transcript_id = row['TranscriptID']
    exon_id = row['ExonID']

    if end_exon and (not allow_incomplete_cds) and end_phase in {1, 2}:
        logging.warn('Exon end phase is 1 or 2 in last exon (%s, %s).',
                     transcript_id, exon_id)
        return False

    if (not end_exon) and end_phase == -1:
        logging.warn('Exon end phase -1 is not in the last exon (%s, %s).',
                     transcript_id, exon_id)
        return False

    return True


def _is_incomplete_cds(row, start_exon, end_exon):
    """
    Return True if there are signals of an incomplete CDS coming from the exon.

    In particular, if the start or end phase is different from 0 or -1.
    """
    if end_exon:
        return row['EndPhase'] in {1, 2}
    if start_exon:
        return row['StartPhase'] in {1, 2}
    return False


def add_protein_seq(data_frame,
                    allow_incomplete_cds=True,
                    seq_column="ExonProteinSequence"):
    """
    Add a column with the protein sequence of the exon.

    It takes a data frame with the exon data, including the sequences as
    the output of add_exon_sequences when it is applied to the output
    of read_exon_file.
    """
    _check_column_presence(data_frame, 'ExonSequence', 'You need to run '
                           'add_exon_sequences first.')
    _check_column_absence(data_frame, seq_column,
                          'Values are going to change.')
    _check_column_absence(data_frame, 'IncompleteCDS',
                          'Values are going to change.')

    sequences = []
    incomplete_cds = []

    row_list = data_frame.to_dict('records')
    n_rows = len(row_list)
    row_index = 0
    while row_index < n_rows:
        start_exon, end_exon = _is_first_or_last_exon(row_list, row_index)

        # Copy of the actual row to have easier access :
        row = row_list[row_index]

        if (row['StartPhase'] == -1) and (row['EndPhase'] == -1):
            if row['Strand'] == -1:
                i = row['ExonRegionEnd'] - row['GenomicCodingEnd']
                j = row['ExonRegionEnd'] - row['GenomicCodingStart']
                j += 1
            else:
                # TO DO: TEST IT!
                i = row['GenomicCodingStart'] - row['ExonRegionStart']
                j = row['GenomicCodingEnd'] - row['ExonRegionStart']
                j += 1

            cds_seq = row['ExonSequence'][i:j]

        else:
            cds_seq = _get_exon_cds(row['ExonSequence'],
                                    row['cDNA_CodingStart'],
                                    row['cDNA_CodingEnd'], start_exon,
                                    end_exon)

            # Ask for the start and end phases of the exon,
            # if the phases are different from 0, the shared codons are
            # assigned to the exon at the left for the negative strand and
            # the exon at the right for the positive strand (i.e. to the
            # exon that changes because of the phase change):

            if not _check_phases_by_position(row, end_exon,
                                             allow_incomplete_cds):
                cds_seq = None

            if row['Strand'] == -1:
                cds_seq = _manage_start_phase_negative_strand(
                    cds_seq, row['StartPhase'])
                cds_seq = _manage_end_phase_negative_strand(
                    row_list, row_index, cds_seq, end_exon)
            else:
                cds_seq = _manage_end_phase_positive_strand(
                    cds_seq, row['EndPhase'])
                cds_seq = _manage_start_phase_positive_strand(
                    row_list, row_index, cds_seq, start_exon)

        # Look for signals of incomplete CDS :
        incomplete = True if cds_seq is None else _is_incomplete_cds(
            row, start_exon, end_exon)

        if not incomplete:
            if (len(cds_seq) % 3) == 0:
                # Add the translated CDS :
                if isinstance(cds_seq, str):
                    sequences.append(
                        Seq(cds_seq, IUPAC.extended_dna).translate())
                else:
                    sequences.append(cds_seq.translate())
            else:
                incomplete = True

        if incomplete:
            sequences.append("")

        incomplete_cds.append(incomplete)

        row_index += 1

    # Add new column to store the sequences :
    data_frame[seq_column] = sequences
    # New column to store signals of incomplete CDS :
    data_frame['IncompleteCDS'] = incomplete_cds


def _different_phases(row_i, row_j):
    """Return True if exon start or end phases are different."""
    value = (row_i['StartPhase'] != row_j['StartPhase']) or \
            (row_i['EndPhase'] != row_j['EndPhase'])
    return value


def _different_coordinates(row_i, row_j):
    """
    Return True if the coordinates values are different.

    Genomic coordinates values are only compared when end and start phases
    are not -1. Phase -1 can be present because of UTR regions, which could
    have different genomic coordinates, but still have the same protein
    sequence.
    """
    if row_i['StartPhase'] != -1 and row_j['StartPhase'] != -1:
        if row_i['GenomicCodingStart'] != row_j['GenomicCodingStart']:
            return True

    if row_i['EndPhase'] != -1 and row_j['EndPhase'] != -1:
        if row_i['GenomicCodingEnd'] != row_j['GenomicCodingEnd']:
            return True

    return False


def _get_comparable_seq(seq_i, seq_j, start_phase, end_phase):
    """
    Return the residues of the exon/subexon that can be shared.

    It takes into account that first/last residue could be different
    if the start/end phase is different from 0 or -1.

    >>> _get_comparable_seq("ASMGSLTSSPSSL", "TSMGSLTSSPSSC", 1, 2)
    ('SMGSLTSSPSS', 'SMGSLTSSPSS')
    """

    if start_phase in [1, 2]:
        seq_i = seq_i[1:]
        seq_j = seq_j[1:]

    if end_phase in [1, 2]:
        seq_i = seq_i[:-1]
        seq_j = seq_j[:-1]

    return seq_i, seq_j


def _identical_seqs(rowi, rowj, seq_column='ExonProteinSequence'):  # pylint: disable=too-many-return-statements
    """Check if the sequences in rows i (rowi) and j (rowj) are identical."""
    if _different_phases(rowi, rowj):
        return False

    if _different_coordinates(rowi, rowj):
        return False

    seq_i = rowi[seq_column]
    seq_j = rowj[seq_column]

    if len(seq_i) != len(seq_j):
        return False

    seq_i, seq_j = _get_comparable_seq(seq_i, seq_j, rowi['StartPhase'],
                                       rowi['EndPhase'])

    return seq_i == seq_j


def find_identical_exons(data_frame,
                         exon_id_column='ExonID',
                         seq_column='ExonProteinSequence'):
    """Find exons that have similar coordinates and identical sequences."""
    exon_clusters = []
    for _, subdf in data_frame.groupby(['GeneID']):
        subdf = subdf.drop_duplicates([
            exon_id_column, 'StartPhase', 'EndPhase', 'GenomicCodingStart',
            'GenomicCodingEnd'
        ])
        row_list = subdf.to_dict('records')
        n_rows = len(row_list)
        for i in range(0, n_rows - 1):
            row_i = row_list[i]
            for j in range(i + 1, n_rows):
                row_j = row_list[j]
                if _identical_seqs(row_i, row_j, seq_column=seq_column):
                    clusters.fill_clusters(exon_clusters,
                                           row_i[exon_id_column],
                                           row_j[exon_id_column])
    return exon_clusters


def _without_stopcodon(sequences):
    """
    Return True for the sequences that doesn't end with '*'.
    """
    sequence = "".join([str(seq) for seq in sequences])
    if sequence:
        return sequence[-1] != "*"
    return True


def delete_incomplete_sequences(data_frame):
    """
    Delete incomplete sequences in place.

    Incomplete sequences are the ones that do not have their protein
    sequence finishing with '*'.
    It also deletes sequences that probably has an incomplete CDS because
    starts or ends with phases differente fron 0 or -1.
    """
    _check_column_presence(data_frame, 'ExonProteinSequence',
                           'You need to run add_protein_seq first.')

    incomplete_cdss = data_frame.loc[:, ['TranscriptID', 'IncompleteCDS'
                                         ]].groupby('TranscriptID').agg(
                                             lambda df: sum(df) > 0)

    incomplete_seqs = data_frame.loc[:, ['TranscriptID',
                                         'ExonProteinSequence']]. \
        groupby('TranscriptID'). \
        agg(_without_stopcodon)  # Incomplete sequences do not end with '*'

    if not incomplete_seqs.empty:
        incomplete_transcripts = set(incomplete_seqs.index[
            incomplete_seqs['ExonProteinSequence']]).union(
                incomplete_cdss.index[incomplete_cdss['IncompleteCDS']])
        if incomplete_transcripts:
            clusters.inform_about_deletions(
                incomplete_transcripts, "Incomplete transcripts were found:")
            data_frame.drop([
                i for i in data_frame.index
                if data_frame.loc[i, 'TranscriptID'] in incomplete_transcripts
            ],
                            inplace=True)


def delete_badquality_sequences(data_frame):
    """Delete protein sequences with X's in their sequence in place."""
    _check_column_presence(data_frame, 'ExonProteinSequence',
                           'You need to run add_protein_seq first.')

    badquality_seqs = data_frame.loc[:, ['TranscriptID',
                                         'ExonProteinSequence']]. \
        groupby('TranscriptID'). \
        agg(lambda df: any('X' in str(s) for s in df))

    if not badquality_seqs.empty:
        badquality_transcripts = set(
            badquality_seqs.index[badquality_seqs['ExonProteinSequence']])
        if badquality_transcripts:
            clusters.inform_about_deletions(
                badquality_transcripts,
                "Transcripts with bad quality sequences were found:")
            data_frame.drop([
                i for i in data_frame.index
                if data_frame.loc[i, 'TranscriptID'] in badquality_transcripts
            ],
                            inplace=True)


def find_identical_sequences(data_frame):
    """
    Find different transcripts with identical sequences.

    Input should have the exon sequences as the output of add_exon_sequences.
    It returns a list with the identical sequence clusters (sets).
    """
    # Concatenate the sequences as a string for each transcript in each gene :
    seqs_df = data_frame.loc[:, ['GeneID',
                                 'TranscriptID',
                                 'ExonProteinSequence']]. \
        groupby(['GeneID', 'TranscriptID']). \
        agg(lambda df: "".join([str(s) for s in df]))
    # Store the cluster of identical sequences :
    seq_clusters = []
    # Find duplicated sequences in each gene transcript :
    for _, subdf in seqs_df.groupby('GeneID'):
        n_rows = len(subdf)
        if n_rows > 1:
            # Sort by 'TranscriptID' to ensure reproducibility
            subdf = subdf.sort_index(level=1)
            row_list = subdf.to_dict('records')
            for i in range(0, n_rows - 1):
                row_i = row_list[i]
                for j in range(i + 1, n_rows):
                    row_j = row_list[j]
                    if row_i['ExonProteinSequence'] == row_j[
                            'ExonProteinSequence']:
                        clusters.fill_clusters(seq_clusters, subdf.index[i][1],
                                               subdf.index[j][1])

    return seq_clusters


def store_cluster(  # pylint: disable=too-many-arguments
        data_frame,
        cluster_list,
        default_values,
        column_name,
        item2str=str,
        get_item=lambda row, default_idx: row[default_idx],
        delim='/'):
    """
    Store the string representation of clusters in the data_frame.

    :param data_frame: DataFrame that's going to store the cluster column.
    :param cluster_list: List of sets.
    :param default_values: Name of existing DataFrame column that containts the default
                           (string) values for the new column.
    :param column_name: Name of the new column to add with the cluster information.
    :param item2str: Function to map from cluster item to its representation in the
                     DataFrame (default: str)
    :param getitem: Function to get the cluster item from the DataFrame row and the index
                    of the default_values column (default: default value).
    :param delim: Character delimiting cluster items (default: '/').


    >>> import pandas as pd
    >>> cluster_list = [{2, 3}]
    >>> df = pd.DataFrame(data={'id': [1, 2, 3, 4], 'value': [1, 2, 2, 3]})
    >>> store_cluster(df, cluster_list, 'id', 'cluster')
    >>> df
       id  value cluster
    0   1      1       1
    1   2      2     2/3
    2   3      2     2/3
    3   4      3       4
    """
    _check_column_absence(data_frame, column_name)
    _check_column_presence(data_frame, default_values)

    cluster_column = []
    default_idx = data_frame.columns.get_loc(default_values) + 1
    # + 1 : the first element of the (row) tuple is the index
    for row in data_frame.itertuples():
        assigned = False
        for cluster in cluster_list:
            if get_item(row, default_idx) in cluster:
                cluster_column.append(
                    clusters.cluster2str(cluster,
                                         delim=delim,
                                         item2str=item2str))
                assigned = True
                break

        if not assigned:
            cluster_column.append(row[default_idx])

    return data_frame.insert(2, column_name, cluster_column)


def delete_identical_sequences(data_frame, seq_column='ExonProteinSequence'):
    """Delete identical sequences in place keeping only one."""
    _check_column_presence(data_frame, seq_column,
                           'You need to run add_protein_seq first.')

    identical_sequences = find_identical_sequences(data_frame)
    store_cluster(data_frame, identical_sequences, 'TranscriptID',
                  'TranscriptIDCluster')
    to_delete = clusters.set_to_delete(identical_sequences)
    clusters.inform_about_deletions(
        to_delete,
        "Identical isoform sequences were found, keeping only one transcript:")
    data_frame.drop([
        i for i in data_frame.index
        if data_frame.loc[i, 'TranscriptID'] in to_delete
    ],
                    inplace=True)


def merge_identical_exons(data_frame,
                          exon_id_column='ExonID',
                          seq_column='ExonProteinSequence'):
    """Unify the 'ExonID' of identical exons."""
    _check_column_presence(data_frame, seq_column,
                           'You need to run add_protein_seq first.')

    identical_exons = find_identical_exons(data_frame,
                                           exon_id_column,
                                           seq_column=seq_column)
    store_cluster(data_frame, identical_exons, exon_id_column,
                  exon_id_column + 'Cluster')
    old2new = dict((exon_id, sorted(group)[0]) for group in identical_exons
                   for exon_id in sorted(group)[1:])
    data_frame.replace(to_replace={exon_id_column: old2new}, inplace=True)


def read_transcript_info(  # pylint: disable=too-many-arguments
        tsl_table_file,
        exon_table_file,
        exon_sequence_file,
        max_tsl_level=3.0,
        remove_na=True,
        remove_badquality=True,
        species_list=None):
    """
    Read and integrate the transcript information.

    Due to the data structure downloaded from ENSEMBL, we combine 3 types of
    primary information:

    1. tsl_table_file has the evidence for the gene transcripts. For each
       transcript, it gives information on its type and its evidence level.
       We use this information to filter out transcripts which are not
       interesting for our analysis (e.g. ncRNA, partial transcripts, etc.)
       At this moment, we keep transcripts with 'Protein coding' biotype and
       with a TSL evidence level of 1, 2 or 3 by default (max_tsl_level).
       If remove_na is True (default) rows with missing TSL data are
       eliminated.
    2. Each row of the exon_table_file summarizes t :whe information for a
       transcript exon, as provided by Biomart. The sorting of the exons
       is important and it is done using the rank information provided by
       ENSEMBL (not coordinates).
    3. The genomic sequence information from the coding sequence of all
       the exons in fasta format is in exon_sequence_file.

    This function also processes the information within transcript to get the
    peptidic sequences of each exon taking into account the exon phases.
    The codons overlapping the exon boundary are always added to the end of
    the first exon. This is different from the older implementation, where it
    is added to the beginning of the next exon.

    This protein sequence information is used by this function to delete
    incomplete and identical sequences and to merge identical exons.

    If remove_badquality is True (default), protein sequences with X's are
    deleted.
    """
    # 1. Read the Transcript Support Level (TSL) table as a pandas' DataFrame:
    tsl_table = read_tsl_file(tsl_table_file,
                              max_tsl_level,
                              remove_na,
                              species_list=species_list)
    # 2. Read the exon table as a pandas' DataFrame:
    exon_table = read_exon_file(exon_table_file)
    # Merge both tables to have all the transcript information in one
    # DataFrame and to filter the exon table based on tsl evidence:
    transcript_info = pd.merge(tsl_table,
                               exon_table,
                               how='inner',
                               on='TranscriptID')
    # Just in case, ensure the correct order after merge:
    transcript_info.sort_values(by=['GeneID', 'TranscriptID', 'ExonRank'],
                                inplace=True)
    # Delete repeated exons:
    transcript_info.drop_duplicates(
        subset=['GeneID', 'TranscriptID', 'ExonID'], inplace=True)
    # 3. Read sequences using BioPython and add them to the exon DataFrame:
    add_exon_sequences(transcript_info, exon_sequence_file)
    # Add protein sequence to the table:
    add_protein_seq(transcript_info)
    # Clean up:
    delete_incomplete_sequences(transcript_info)
    if remove_badquality:
        delete_badquality_sequences(transcript_info)

    if transcript_info.empty:
        raise ValueError("There are no transcripts left!")

    # Ensure exon order and check the downloaded exon phase information:
    phases.check_order_and_phases(transcript_info)
    # NOTE: delete_incomplete_sequences deletes sequences with incomplete CDSs.
    # check_order_and_phases should run after it, because it throws errors with
    # some incomplete CDSs.
    # Delete incomplete CDS identified by check_order_and_phases:
    delete_incomplete_sequences(transcript_info)
    # More clean up:
    delete_identical_sequences(transcript_info)
    merge_identical_exons(transcript_info)

    return transcript_info
