"""
phases: Calculate exon phases and chack exon order and phases.

This is the code for calculating exon start and end phases that we used before
downloading that data from ENSEMBL. The code was changed to be use as a checker
for the downloaded and parsed data. Also, this module have functions used to
calculate the subexon phases.

This could be useful to help in the understanding of exon start and end phases.

Static function call graph:
---------------------------
::

    +----------------------------------+
    |      check_order_and_phases      |
    +-+--------+-----------------------+
      |        |
      |        |
      |        v
      |      +-+-----------------------+
      |      |    _check_exon_order    |
      |      +-------------------------+
      v
    +-+--------------------------------+
    |          _check_phases           |
    +-+--------+---------+-------------+
      |        |         |
      |        |         |
      |        |         v
      |        |       +-+-------------+
      |        |       | _equal_phases |
      |        v       +---------------+
      |      +-+-----------------------+
      |      | end_phase_previous_exon |
      ^      +-------------------------+
    +-+--------------------------------+
    |         calculate_phase          |
    +-+--------------------------------+
      |
      |
      v
    +-+--------------------------------+
    | bases_to_complete_previous_codon |
    +----------------------------------+

"""

import logging


def bases_to_complete_previous_codon(phase):  # pylint: disable=invalid-name
    """
    Return the bases at exon start that completes the previous codon.

    >>> bases_to_complete_previous_codon(0) # e.g. -----XXX
    0
    >>> bases_to_complete_previous_codon(2) # e.g. XX-----X
    1
    >>> bases_to_complete_previous_codon(1) # e.g. X-----XX
    2
    """
    if phase == 0:
        return 0
    if phase == 1:
        return 2
    if phase == 2:
        return 1

    logging.warning(
        "Only phases 0, 1 and 2 are allowed for internal exons (phase: %s)" %
        phase)
    return None


def bases_to_complete_next_codon(phase):  # pylint: disable=invalid-name
    """
    Return the bases at the exon end that completes the next codon.

    >>> bases_to_complete_next_codon(0) # e.g. -----XXX
    0
    >>> bases_to_complete_next_codon(2) # e.g. XX-----X
    2
    >>> bases_to_complete_next_codon(1) # e.g. X-----XX
    1
    """
    if phase in {0, 1, 2}:
        return phase

    logging.warning(
        "Only phases 0, 1 and 2 are allowed for internal exons (phase: %s)" %
        phase)
    return None


def calculate_phase(cdna_len, previous_end_phase):
    """
    Calculate exon start and end phases using the intron phase.

    It calculates phases as explained in the ENSEMBL glossary:
    http://www.ensembl.org/Multi/Help/Glossary

    "In protein-coding exons, the end phase is the place where the intron
    lands inside the codon : 0 between codons, 1 between the 1st and second
    base, 2 between the second and 3rd base. Exons therefore have a start
    phase and an end phase, but introns have just one phase. An exon which
    is non coding (or non-coding at the end) has an end phase of -1"

    Under this definition, start phase of the actual exon is the end phase
    of the previous one, and the end phase of the actual is the start phase
    of the next one.

    Let X, y and z denote codon bases that belong to exons X, y and z.
    Let - denote intron bases.

    0 - No interruption.
        XXX------yyyzzz

    1 - first codon's first base is in the previous exon.
        X------XXyyyzzz

    2 - first codon's first two bases are in the
        previous exon.
        XX------Xyyyzzz

    NOTE: This function can not calculate a phase of -1.

    >>> calculate_phase(6, 0) # e.g. -----XXXyyy----zzz
    (0, 0)
    >>> calculate_phase(6, 2) # e.g. XX-----Xyyyzz----z
    (2, 2)
    >>> calculate_phase(6, 1) # e.g. X-----XXyyyz----zz
    (1, 1)
    """
    n_bases = bases_to_complete_previous_codon(previous_end_phase)
    if n_bases is None:
        return None, None
    next_start_phase = (
        cdna_len - n_bases) % 3

    return previous_end_phase, next_start_phase


def end_phase_previous_exon(data_frame,
                            exon_pos,
                            prev_row_index,
                            end_phase_column='EndPhase'):
    """
    Return the end phase of the previous exon.

    It returns 0 if the actual exon is the first in the transcript.
    """
    if exon_pos == 0:
        return 0

    return data_frame.loc[prev_row_index, end_phase_column]


def _equal_phases(downloaded_phase, calculated_phase):
    """
    Compare calculated and downloaded phases.

    If the downloaded_phase is -1, this function returns True, otherwise,
    the phases are compared and True is returned if they are equal.

    >>> _equal_phases(0, 0)
    True
    >>> _equal_phases(-1, 1)
    True
    >>> _equal_phases(2, 1)
    False
    """
    if downloaded_phase != -1:
        if downloaded_phase != calculated_phase:
            return False

    return True


def _check_phases(data_frame, row_number, row_index, prev_row_index, exon_pos):
    """
    Check the calculated phases against the downloaded ones.

    It takes into account that -1 phases cannot be calculated, and therefore
    they are not compared (thanks to _equal_phases).

    Return True if everything is OK.
    """
    exon_id = data_frame.loc[row_index, 'ExonID']
    transcript_id = data_frame.loc[row_index, 'TranscriptID']

    previous_end_phase = end_phase_previous_exon(data_frame, exon_pos,
                                                 prev_row_index)

    cdna_len = data_frame.loc[row_index, 'cDNA_CodingEnd'] - \
        data_frame.loc[row_index, 'cDNA_CodingStart'] + 1

    start_phase, end_phase = calculate_phase(cdna_len, previous_end_phase)
    if (start_phase is None) or (end_phase is None):
        return False

    df_start_phase = data_frame.loc[row_index, 'StartPhase']
    df_end_phase = data_frame.loc[row_index, 'EndPhase']

    if not _equal_phases(df_start_phase, start_phase):
        logging.warning(
            "%d start phase is observed but %d is expected (%s, %s).",
            df_start_phase, start_phase, transcript_id, exon_id)
        return False

    if (df_start_phase == -1) and (exon_pos != 0):
        logging.warning("start phase -1 is not in the first exon (%s, %s).",
                        transcript_id, exon_id)
        return False

    if not _equal_phases(df_end_phase, end_phase):
        logging.warning(
            "%d start phase is observed but %d is expected (%s, %s).",
            df_end_phase, end_phase, transcript_id, exon_id)
        return False

    n_rows = data_frame.shape[0]
    if df_end_phase == -1 and row_number != n_rows - 1:
        # next row index.
        next_row_index = data_frame.index[row_number + 1]
        if transcript_id == data_frame.loc[next_row_index, 'TranscriptID']:
            logging.warning("end phase -1 is not in the last exon (%s, %s).",
                            transcript_id, exon_id)
            return False

    return True


def _check_exon_order(data_frame, row_index, prev_row_index, exon_pos):
    """
    Check that exons are ordered by rank in the transcript.
    
    Return True if everything is OK.
    """
    exon_id = data_frame.loc[row_index, 'ExonID']
    transcript_id = data_frame.loc[row_index, 'TranscriptID']
    if (exon_pos != 0) and (data_frame.loc[row_index, 'ExonRank'] <=
                            data_frame.loc[prev_row_index, 'ExonRank']):
        logging.warning("Exons aren't sorted by rank (%s, %s).", transcript_id,
                        exon_id)
        return False
    return True


def check_order_and_phases(data_frame):
    """
    Check DataFrame order and exon start/end phases.

    It takes a data frame with the exon data, including the sequences as
    the output of add_exon_sequences when it is applied to the output
    of read_exon_file.

    The columns of the input dataframe must be ordered by 'Transcript
    stable ID' and 'ExonRank'. The function read_exon_file
    ensures that order.

    It checks that order and it also checks exon start and end phase
    information. It throws an informative error if something is wrong.
    """
    # Use a set of transcript ids to ensure that data is sorted
    finished_transcripts = set()

    assert data_frame.loc[data_frame.index[0], 'StartPhase'] in {-1, 0},\
        'First exon has start phase of 1 or 2, row: %s' % data_frame.iloc[0, :]

    n_rows = data_frame.shape[0]
    row_number = 1  # to skip the first row
    exon_pos = 1  # to store exon rank without counting non coding exons
    # as 'ExonRank' does.
    while row_number < n_rows:

        # data_frame.loc[index, column] uses the index of the row (i)
        # rather than the integer location (row_number) :
        row_index = data_frame.index[row_number]
        # previous row index.
        prev_row_index = data_frame.index[row_number - 1]
        # If transcript id changes, the next exon is 0,
        # i.e. the first of the next gene transcript :
        actual_transcript = data_frame.loc[row_index, 'TranscriptID']
        prev_transcript = data_frame.loc[prev_row_index, 'TranscriptID']

        if actual_transcript != prev_transcript:
            finished_transcripts.add(prev_transcript)
            exon_pos = 0
        else:
            assert actual_transcript not in finished_transcripts, \
                "Data isn't sorted by 'TranscriptID'"
            exon_pos += 1

        if (not _check_exon_order(data_frame, row_index, prev_row_index,
                                  exon_pos)) or (not _check_phases(
                                      data_frame, row_number, row_index,
                                      prev_row_index, exon_pos)):
            data_frame.at[row_index, 'IncompleteCDS'] = True

        row_number += 1
