import pandas as pd
import numpy as np
import pytest

from thoraxe.subexons.tidy import get_tidy_table


def _mk_row(gene="G", trx="T", subexon_rank=15, strand=-1, start=0, end=None, label="X"):
    # End defaults to one base before start on the negative strand or one base after on positive
    if end is None:
        end = start - 3 if strand == -1 else start + 3

    # Keep protein length consistent with S_exon_Lengths ("21")
    subexon_prot = "A" * 21
    subexon_nt = "ATG" * 7  # 21 aa -> 63 nt; not strictly required but tidy to keep

    # NOTE: get_tidy_table reads S_exon_Sequences (plural) from input rows
    s_exon_seq = "A" * 21

    return {
        "GeneID": gene,
        "TranscriptIDCluster": trx,
        "SubexonRank": subexon_rank,
        "Strand": strand,
        "SubexonCodingStart": int(start),
        "SubexonCodingEnd": int(end),

        # Columns consumed by tidy.get_tidy_table()
        "SubexonProteinSequence": subexon_prot,
        "SubexonSequence": subexon_nt,

        # One s-exon per row in these tests; order is what we care about
        "S_exons": label,
        "S_exon_Lengths": "21",
        "S_exon_Sequences": s_exon_seq,
    }


@pytest.mark.parametrize(
    "species,trx,strand,starts",
    [
        # Platypus (ornithorhynchus_anatinus) transcripts
        ("ornithorhynchus_anatinus", "ENSOANT00000050045", -1,
         {"17_2": 29149002, "1_0": 29143989, "1_1": 29143910}),
        ("ornithorhynchus_anatinus", "ENSOANT00000062336", -1,
         {"17_2": 29149002, "1_0": 29143989, "1_1": 29143910}),
        # Xenopus (xenopus_tropicalis) transcript
        ("xenopus_tropicalis", "ENSXETT00000106354", -1,
         {"17_2": 3169580, "1_0": 3164877, "1_1": 3164798}),
    ],
)
def test_real_negative_strand_order(species, trx, strand, starts):
    # Create three rows with the same SubexonRank to trigger tie-breaking
    # Intentionally insert rows in an order that is NOT genomic order
    # so the test fails if get_tidy_table does not break ties using
    # genomic start positions (via the signed-start sort key).
    insertion_order = ["1_1", "1_0", "17_2"]
    rows = [
        _mk_row(
            gene="G",
            trx=trx,
            subexon_rank=15,
            strand=strand,
            start=starts[label],
            end=starts[label] - 1,
            label=label,
        )
        for label in insertion_order
    ]
    pre = pd.DataFrame(rows)

    # gene2species mapping is required by get_tidy_table in your ThorAxe branch
    tidy = get_tidy_table(pre, gene2species={"G": species})

    grp = tidy[tidy["TranscriptIDCluster"] == trx].sort_values("S_exon_Rank")
    got = grp["S_exonID"].tolist()

    # Expected order: for negative strand, 5′→3′ is decreasing start.
    # The tidy table exposes `SubexonCodingStart` (not `S_exon_CodingStart`).
    signed = np.where(
        grp["Strand"].to_numpy() == 1,
        grp["SubexonCodingStart"].to_numpy(),
        -grp["SubexonCodingStart"].to_numpy(),
    )
    exp = grp.iloc[np.argsort(signed, kind="stable")]["S_exonID"].tolist()

    assert got == exp, f"{species}:{trx} expected {exp} by genomic order, got {got}"
