"""
tidy:  Utils for having a tidy output (denormalized table).
"""

import os

import pandas as pd

from itertools import chain
from thoraxe.utils import split


def get_tidy_table(table, gene2species):
    """
    Save a tidy and denormalized output table.

    Takes the final exon table and a dict from gene to species name to save a
    csv file with the output table in the outputdir.
    """
    table = table.sort_values(
        by=["GeneID", "TranscriptIDCluster", "SubexonRank"])

    tidy_table = pd.DataFrame(columns=[
        'Species', "GeneID", "TranscriptIDCluster", "ExonIDCluster",
        "SubexonIDCluster", "HomologousExon", "Strand", "TranscriptsInGene",
        "ExonRank", "Cluster", "StartPhase", "EndPhase", "SubexonRank",
        "SubexonSequence", "SubexonProteinSequence", "SubexonCodingStart",
        "SubexonCodingEnd", "TranscriptsWithSubexon", "HomologousExonRank",
        "HomologousExonSequence", "HomologousExonStart", "HomologousExonEnd"
    ])

    previous_transcript = ""
    for row in table.itertuples():
        if previous_transcript != row.TranscriptIDCluster:
            previous_transcript = row.TranscriptIDCluster
            homologous_exon_rank = 0

        sube_seq = str(row.SubexonProteinSequence).replace('*', '')
        hexons = split.split_exons(row.HomologousExons)
        he_lens = split.split_lengths(row.HomologousExonLengths, sube_seq)
        he_seqs = split.split_seqs(row.HomologousExonSequences, sube_seq)
        assert len(hexons) == len(he_lens)
        assert len(hexons) == len(he_seqs)

        tidy_row = dict()
        for column in tidy_table.columns:
            if column in table.columns:
                tidy_row[column] = getattr(row, column)

        tidy_row['Species'] = gene2species[row.GeneID]

        stop = 0
        for (hexon, he_len, he_seq) in zip(hexons, he_lens, he_seqs):
            homologous_exon_rank += 1
            tidy_row['Species'] = gene2species[row.GeneID]
            tidy_row['HomologousExon'] = hexon
            tidy_row['HomologousExonRank'] = homologous_exon_rank
            tidy_row['HomologousExonSequence'] = he_seq
            if stop == 0:
                start = 1
                stop = he_len
            else:
                start = stop + 1
                stop = stop + he_len
            tidy_row['HomologousExonStart'] = start
            tidy_row['HomologousExonEnd'] = stop
            tidy_table = tidy_table.append(tidy_row, ignore_index=True)

    return tidy_table
