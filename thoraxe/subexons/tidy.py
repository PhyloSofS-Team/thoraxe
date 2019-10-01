"""
tidy:  Utils for having a tidy CSV table output (denormalized table).
"""

import pandas as pd

from thoraxe.utils import split


def get_tidy_table(table, gene2species):  # pylint: disable=too-many-locals
    """
    Save a tidy and denormalized output table.

    Takes the final exon table and a dict from gene to species name to save a
    csv file with the output table in the outputdir.
    """
    table = table.sort_values(
        by=["GeneID", "TranscriptIDCluster", "SubexonRank"])

    tidy_table = pd.DataFrame(columns=[
        'Species', "GeneID", "TranscriptIDCluster", "ExonIDCluster",
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
        he_lens = split.split_lengths(row.S_exon_Lengths, sube_seq)
        he_seqs = split.split_seqs(row.S_exon_Sequences, sube_seq)
        assert len(hexons) == len(he_lens)
        assert len(hexons) == len(he_seqs)

        tidy_row = dict()
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
            tidy_table = tidy_table.append(tidy_row, ignore_index=True)

    return tidy_table
