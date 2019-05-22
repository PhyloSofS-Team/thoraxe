"""
tidy:  Utils for having a tidy output (denormalized table).
"""

import os

from itertools import chain
from thoraxe.utils import split


def save_tidy_csv(outputdir, data, gene2species):
    """
    Save a tidy and denormalized output table.

    Takes the final exon table and a dict from gene to species name to save a
    csv file with the output data in the outputdir.
    """
    first_columns = ['Species']
    columns = [
        "GeneID", "TranscriptIDCluster", "ExonIDCluster", "SubexonIDCluster",
        "HomologousExon", "Strand", "TranscriptsInGene", "ExonRank",
        "ExonProteinSequence", "Cluster", "StartPhase", "EndPhase",
        "SubexonRank", "SubexonSequence", "SubexonCodingStart",
        "SubexonCodingEnd", "TranscriptsWithSubexon"
    ]
    last_columns = [
        "HomologousExonRank",
        "HomologousExonSequence",
        "HomologousExonStart",
        "HomologousExonEnd",
    ]

    data = data.sort_values(
        by=["GeneID", "TranscriptIDCluster", "SubexonRank"])

    with open(os.path.join(outputdir, "homologous_exon_table.csv"),
              "w") as outfile:
        for column in chain(first_columns, columns, last_columns[:-1]):
            outfile.write(column)
            outfile.write(',')
        outfile.write(last_columns[-1])
        outfile.write('\n')

        previous_transcript = ""
        for row in data.itertuples():
            if previous_transcript != row.TranscriptIDCluster:
                previous_transcript = row.TranscriptIDCluster
                homologous_exon_rank = 0

            sube_seq = str(row.SubexonSequence).replace('*', '')
            hexons = split.split_exons(row.HomologousExons)
            he_lens = split.split_lengths(row.HomologousExonLengths, sube_seq)
            he_seqs = split.split_lengths(row.HomologousExonSequences,
                                          sube_seq)
            assert len(hexons) == len(he_lens)
            assert len(hexons) == len(he_seqs)
            outfile.write(gene2species[row.GeneID])
            outfile.write(',')
            for column in columns:
                outfile.write(getattr(row, column))
                outfile.write(',')

            stop = 0
            for (hexon, he_len, he_seq) in zip(hexons, he_lens, he_seqs):
                homologous_exon_rank += 1
                outfile.write(homologous_exon_rank)
                outfile.write(',')
                outfile.write(he_seq)
                outfile.write(',')
                if stop == 0:
                    start = 1
                    stop = he_len
                else:
                    start = stop + 1
                    stop = stop + he_len
                outfile.write(start)
                outfile.write(',')
                outfile.write(stop)
                outfile.write('\n')
