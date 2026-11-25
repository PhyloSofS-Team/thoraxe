"""
tidy:  Utils for having a tidy CSV table output (denormalized table).
"""

import pandas as pd

from thoraxe.utils import split


def _own_s_exon_sequences(lengths, own):
    """
    Cut the row's own sub-exon protein at the s-exon lengths.

    The s-exon sequences stored per sub-exon come from one representative
    row, so a codon split at a sub-exon boundary can carry another
    transcript's residue, or one residue more or less. Each row gets its own
    residues instead; a length difference goes to the first s-exon (the
    sub-exon start), or to the last one if the first is too short.

    >>> _own_s_exon_sequences([2, 3], 'MKLAV')
    ([2, 3], ['MK', 'LAV'])
    >>> _own_s_exon_sequences([2, 3], 'MKLA')
    ([1, 3], ['M', 'KLA'])
    """
    lengths = [int(length) for length in lengths]
    diff = len(own) - sum(lengths)
    if diff:
        index = 0 if lengths[0] + diff >= 0 else len(lengths) - 1
        lengths[index] = max(0, lengths[index] + diff)
    seqs, pos = [], 0
    for length in lengths:
        seqs.append(own[pos:pos + length])
        pos += length
    return lengths, seqs


def get_tidy_table(table, gene2species):  # pylint: disable=too-many-locals
    """
    Save a tidy and denormalized output table.

    Takes the final exon table and a dict from gene to species name to save a
    csv file with the output table in the outputdir.
    """
    # Add a signed start column to sort by start position taking into account
    # the strand
    table = table.assign(
        _SignedStart=table["SubexonCodingStart"] * table["Strand"]
    )
    # Sort by GeneID, TranscriptIDCluster, SubexonRank and signed start
    # (to have a consistent order of s-exons when they have the same rank)
    table = table.sort_values(
        by=["GeneID", "TranscriptIDCluster", "SubexonRank", "_SignedStart"],
        kind="stable"
    )
    # Remove the signed start column
    table = table.drop(columns=["_SignedStart"])

    tidy_table = pd.DataFrame(columns=[
        "Species", "GeneID", "TranscriptIDCluster", "ExonIDCluster",
        "SubexonIDCluster", "S_exon", "Strand", "TranscriptsInGene",
        "ExonRank", "Cluster", "StartPhase", "EndPhase", "SubexonRank",
        "SubexonSequence", "SubexonProteinSequence", "SubexonCodingStart",
        "SubexonCodingEnd", "TranscriptsWithSubexon", "S_exon_Rank",
        "S_exon_Sequence", "S_exon_Start", "S_exon_End"
    ])

    previous_transcript = ""
    for row in table.itertuples():
        if previous_transcript != row.TranscriptIDCluster:
            previous_transcript = row.TranscriptIDCluster
            s_exon_rank = 0
        sube_seq = str(row.SubexonProteinSequence).replace('*', '')
        hexons = split.split_exons(row.S_exons)
        if sube_seq:
            he_lens = split.split_lengths(row.S_exon_Lengths, sube_seq)
            he_seqs = split.split_seqs(row.S_exon_Sequences, sube_seq)
            assert len(hexons) == len(he_lens)
            assert len(hexons) == len(he_seqs)
            he_lens, he_seqs = _own_s_exon_sequences(he_lens, sube_seq)
        else:
            he_lens = [0]
            he_seqs = ['']

        tidy_row = {}
        for column in list(tidy_table.columns):
            if column in list(table.columns):
                tidy_row[column] = getattr(row, column)

        tidy_row['Species'] = gene2species[row.GeneID]

        stop = 0
        for (hexon, he_len, he_seq) in zip(hexons, he_lens, he_seqs):
            s_exon_rank += 1
            tidy_row['Species'] = gene2species[row.GeneID]
            tidy_row['S_exon'] = hexon
            tidy_row['S_exon_Rank'] = s_exon_rank
            tidy_row['S_exon_Sequence'] = he_seq
            if stop == 0:
                start = 1
                stop = he_len
            else:
                start = stop + 1
                stop = stop + he_len
            tidy_row['S_exon_Start'] = start
            tidy_row['S_exon_End'] = stop
            # FIX for warning in future versions (03/2024)
            #tidy_table = tidy_table.append(tidy_row, ignore_index=True)
            tidy_table = pd.concat([tidy_table,pd.DataFrame([tidy_row])], ignore_index=True)

    tidy_table.rename(columns={
        "StartPhase": "SubexonStartPhase",
        "EndPhase": "SubexonEndPhase",
        "S_exon": "S_exonID"
    },
                      inplace=True)

    return tidy_table
