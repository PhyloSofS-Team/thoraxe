import os
import pytest

from thoraxe import transcript_info
from thoraxe import subexons


@pytest.fixture(scope='module')
def clustered_trx_data(request):
    filename = request.module.__file__

    def _get_clustered_trx_data(folder):
        test_dir = os.path.dirname(filename)
        datadir = os.path.join(test_dir, 'data')
        folder_path = os.path.join(datadir, folder, 'Ensembl')
        trx_data = transcript_info.read_transcript_info(
            os.path.join(folder_path, 'tsl.csv'),
            os.path.join(folder_path, 'exonstable.tsv'),
            os.path.join(folder_path, 'sequences.fasta'),
            remove_na=False)
        return transcript_info.exon_clustering(trx_data)

    return _get_clustered_trx_data


def _st(l):
    """Get the first element."""
    return list(l)[0]


def test_subexon_table(clustered_trx_data):
    data = clustered_trx_data('GPRIN1')

    # There are genes in both strands
    assert sum(data.Strand == 1) == 2
    assert sum(data.Strand == -1) == 7
    # GenomicCodingStart of a negative-strand gene is its GenomicCodingEnd at
    # this point. The final sub-exon table correct this:
    assert sum(data.GenomicCodingStart > data.GenomicCodingEnd) == 0

    # negative-strand gene
    human_seq = ("MDTAEDPAWLQLLQKDSSPPGPRPTAFFCPQDGSLGAGSSAMRDYCPSQQKASPAPPRHT"
                 "PDQSPGMESRHRSPSGAGEGASCSDGPRGSLACPSPTCFSPQESPSKETLEAHGASISGT"
                 "PEATTSGKPEPVSSVKTEPKSSDDRNPMFLEKMDFKSSKQADSTSIGKEDPGSSRKADPM"
                 "FTGKAEPEILGKGDPVAPGRMDPMTVRKEDLGSLGKVDPLCSSKTYTVSPRKEDPGSLRK"
                 "VDPVSSDKVDPVFPRKEEPRYSGKEHPVSSEKVAPTSAEKVDLVLSGKRDPGPSGKADPM"
                 "PLESMDSASTGKTEPGLLGKLIPGSSGKNGPVSSGTGAPGSLGRLDPTCLGMADPASVGN"
                 "VETVPATKEDSRFLGKMDPASSGEGRPVSGHTDTTASAKTDLTSLKNVDPMSSGKVDPVS"
                 "LGKMDPMCSGKPELLSPGQAERVSVGKAGTVSPGKEDPVSSRREDPISAGSRKTSSEKVN"
                 "PESSGKTNPVSSGPGDPRSLGTAGPPSAVKAEPATGGKGDPLSSEKAGLVASGKAAPTAS"
                 "GKAEPLAVGKEDPVSKGKADAGPSGQGDSVSIGKVVSTPGKTVPVPSGKVDPVSLGKAEA"
                 "IPEGKVGSLPLEKGSPVTTTKADPRASGKAQPQSGGKAETKLPGQEGAAAPGEAGAVCLK"
                 "KETPQASEKVDPGSCRKAEPLASGKGEPVSLGKADSAPSRKTESPSLGKVVPLSLEKTKP"
                 "SSSSRQLDRKALGSARSPEGARGSEGRVEPKAEPVSSTEASSLGQKDLEAAGAERSPCPE"
                 "AAAPPPGPRTRDNFTKAPSWEASAPPPPREDAGTQAGAQACVSVAVSPMSPQDGAGGSAF"
                 "SFQAAPRAPSPPSRRDAGLQVSLGAAETRSVATGPMTPQAAAPPAFPEVRVRPGSALAAA"
                 "VAPPEPAEPVRDVSWDEKGMTWEVYGAAMEVEVLGMAIQKHLERQIEEHGRQGAPAPPPA"
                 "ARAGPGRSGSVRTAPPDGAAKRPPGLFRALLQSVRRPRCCSRAGPTAE*")

    human_row = data.ExonIDCluster == "ENSE00001159289"
    assert str(_st(data.loc[human_row, "ExonProteinSequence"])) == human_seq
    assert _st(data.loc[human_row, "StartPhase"]) == -1
    assert _st(data.loc[human_row, "EndPhase"]) == -1
    assert _st(data.loc[human_row, "GenomicCodingStart"]) == 176596808
    assert _st(data.loc[human_row, "GenomicCodingEnd"]) == 176599834
    assert _st(data.loc[human_row, "GenomicCodingStart"]) > 176595802
    assert _st(data.loc[human_row, "GenomicCodingEnd"]) < 176599877

    # positive-strand gene
    pig_seq = ("MKERCPSQQRASPAPSRQTTTQSPGMDSRHSSPIGAGEGASCSEGPDGSLAYPSLTCILP"
               "QEAAAKETLGTDGALISGTPETTLSGKPEPVSSVKTGISSSENRNSLFLEKMDSKSSKQA"
               "DAISIGKEDAGSLRKADPMFTGKTEPAILGKGAPVASERMDPMTPRKEDPGSLGKVDPEC"
               "SDKIDTMAPRREDPGSLGKVELAFPREEEPGNSGEELPVSSEEVSPTSVSKTALISFGKR"
               "DPESSGKVGPMTLENMDSASAGKTDPGSLGSTPGSSGKIKPVSPGTKALGSSERVDPTSL"
               "GMTDPVSVLNAETVSSTKEDPQVLGKMDPASSGKQDPMSMRMTKTVSVGQVDPMFSGKMD"
               "TTSLKNKEPVSSGKVSPVSLGKVDPVPSGKPETLSPGKAEPTSGGKTETASSGKEDPVSS"
               "GKVDPSSVGNIKTSSSGKVNPESSGKMDPVSSGPGGPPSSVKAEAVTGGKADLLPLEKSG"
               "PVAPGKVGPVTSGKADSPTLGKVDPVSKGKTETVAPGDKDPVSAGKVAPTTLGGTVPASL"
               "GKAEAVPEGKVNPLPLEKDNPMNSTKLDPGASGKAEPKSEGKAETKPPGQEGTASSGKEE"
               "AKSLQKEKPLTLEKGDPGSSGKADPPASGEPAALGKTDSVPPRKAEPASSEKEAHLTLQK"
               "EEVSSSRQADDKPCGSAPSAAGLGGRRDPAPPPSSQTSSLSQKDPAAARAESSPGLEAVG"
               "QPPGPRTRDNFTKAPSWEASAPPPREDAGTQAGAQACVSVAVSPMSPQDGAGGPAFSFQA"
               "KQRAPSPAPGPPSRRDAGLQVSLGAAETRSVATGPMTPQAAAPPAFPEVRVRPGSALAAA"
               "VAPPEAAEPVRDVSWDEKGMTWEVYGAAMEVEVLGMAIQKHLERQIEEHGRQGAPAPPPA"
               "ARPGLGTGRAGSVRAAPPRAPPNARRAFSAHCCRVCAGRGAARGPDPRPSDQLRFCTPGF"
               "PTFSGSLLDRRPLEGIPSIHDRAARGRAASRVFSSSRLSSPLPS*")

    pig_row = data.ExonIDCluster == "ENSSSCE00000122667"
    assert str(_st(data.loc[pig_row, "ExonProteinSequence"])) == pig_seq
    assert _st(data.loc[pig_row, "StartPhase"]) == 0
    assert _st(data.loc[pig_row, "EndPhase"]) == 0
    assert _st(data.loc[pig_row, "GenomicCodingStart"]) == 81315368
    assert _st(data.loc[pig_row, "GenomicCodingEnd"]) == 81318382

    subexon_table = subexons.create_subexon_table(data)

    pig_se = subexon_table.SubexonID == "ENSSSCE00000122667_SE_0"
    human_se = subexon_table.SubexonID == "ENSE00001159289_SE_0"

    assert str(_st(
        subexon_table.loc[human_se, "SubexonProteinSequence"])) == human_seq
    assert _st(subexon_table.loc[human_se, "Strand"]) == -1
    assert _st(subexon_table.loc[human_se, "StartPhase"]) == 0
    assert _st(subexon_table.loc[human_se, "EndPhase"]) == 0
    # Correct order:
    assert _st(subexon_table.loc[human_se, "SubexonCodingStart"]) == 176599834
    assert _st(subexon_table.loc[human_se, "SubexonCodingEnd"]) == 176596808

    assert str(_st(
        subexon_table.loc[pig_se, "SubexonProteinSequence"])) == pig_seq
    assert _st(subexon_table.loc[pig_se, "Strand"]) == 1
    assert _st(subexon_table.loc[pig_se, "StartPhase"]) == 0
    assert _st(subexon_table.loc[pig_se, "EndPhase"]) == 0
    assert _st(subexon_table.loc[pig_se, "SubexonCodingStart"]) == 81315368
    assert _st(subexon_table.loc[pig_se, "SubexonCodingEnd"]) == 81318382

    # SubexonCodingStart is correct :
    assert sum(
        subexon_table.SubexonCodingStart > subexon_table.SubexonCodingEnd) == 7
    assert sum(
        subexon_table.SubexonCodingStart < subexon_table.SubexonCodingEnd) == 2

    assert 'SubexonID' in subexon_table.columns
    assert 'SubexonIDCluster' in subexon_table.columns

    assert 'ENSMUST00000099506/ENSMUST00000135343' in subexon_table[
        'TranscriptIDCluster'].values  # .values, otherwise in
    # checks whether the value is in the Series' index.

    assert all(subexon_table['TranscriptFraction'] == 1.0)


def test_subexon_clusters(clustered_trx_data):
    data = clustered_trx_data('MAPK8')
    subexon_table = subexons.create_subexon_table(data)
    assert len(
        subexon_table.loc[subexon_table['SubexonIDCluster'] ==
                          'ENSMMUE00000040028_SE_0/ENSMMUE00000387789_SE_0',
                          'SubexonID'].unique()) == 1

    not_merged = subexons.create_subexon_table(data, merge_non_redundant=False)
    assert not_merged.shape[0] > subexon_table.shape[0]

    # QVQQ
    assert len(
        not_merged[not_merged['SubexonID'] == 'ENSMUSE00000689835_SE_2']) > 0
    assert len(
        not_merged[not_merged['SubexonID'] == 'ENSMUSE00000689835_SE_1']) > 0
    assert len(subexon_table[subexon_table['SubexonID'] ==
                             'ENSMUSE00000689835_SE_2']) == 0
    assert len(subexon_table[subexon_table['SubexonID'] ==
                             'ENSMUSE00000689835_SE_1']) == 0
    qvqq = subexon_table[subexon_table['SubexonID'] ==
                         'ENSMUSE00000689835_SE_2_1']
    assert len(qvqq) == 2
    assert sorted(
        qvqq['TranscriptID']) == ['ENSMUST00000111943', 'ENSMUST00000111945']
    for index in [0, 1]:
        assert qvqq['Strand'].iloc[index] == -1
        assert str(qvqq['SubexonProteinSequence'].iloc[index]) == 'QVQQ*'
        assert str(qvqq['SubexonSequence'].iloc[index]) == 'CACAGGTGCAGCAATGA'
        assert qvqq['ExonIDCluster'].iloc[
            index] == 'ENSMUSE00000689835/ENSMUSE00000689841'

    # VINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR
    exon_row = subexon_table[subexon_table['ExonID'] ==
                             'ENSMMUE00000040064'].iloc[0]
    assert exon_row['SubexonID'] == 'ENSMMUE00000040064_SE_15_16'
    assert exon_row['SubexonIDCluster'] == 'ENSMMUE00000040064_SE_15_16'
    assert str(exon_row['SubexonProteinSequence']
               ) == 'AVINGSQHPSSSSSVNDVSSMSTDPTLASDTDSSLEASAGPLGCCR*'
    assert str(exon_row['SubexonSequence']) == (
        'CAGTGATCAATGGCTCTCAGCATCCATCGTCATCGTCGTCTGTCAATGATGTGTCTTCAA'
        'TGTCAACAGATCCGACTTTGGCCTCGGATACAGACAGCAGTCTAGAAGCATCAGCTGGGC'
        'CTCTGGGCTGCTGTAGATGA')
    assert exon_row['Strand'] == 1
    assert exon_row['SubexonCodingStart'] == 42539779
    assert exon_row['SubexonCodingEnd'] == 42539918

    # Coordinates and phases

    # ENSGGOE00000120742 has two sub-exons: se_a, se_b

    # se_b :
    # \/ 62034360 (phase: 1) coding start
    # ATATTGATCAGTGGAATAAAGTTATTGAACAGCTTGGAACACCATGTCCTGAATTCATGA
    # AGAAGCTGCAACCAACAGTAAGGACTTACGTTGAAAACAGACCTAAATATGCTGGATATA
    # GCTTTGAGAAACTCTTCCCTGATGTCCTTTTCCCAG
    #                                    /\ 62034205 (phase: 1)
    # se_a :
    # CTGACTCAGAACACAACAAACTTAAAG
    #                          /\ 62034178 (phase: 1) coding end
    # /\ 62034204 (phase: 1)

    se_a = subexon_table.loc[156, :]
    se_b = subexon_table.loc[169, :]

    assert se_a["SubexonID"] == "ENSGGOE00000120742_SE_7"
    assert se_b["SubexonID"] == "ENSGGOE00000120742_SE_8"

    assert se_a[
        "SubexonIDCluster"] == "ENSGGOE00000120742_SE_7/ENSGGOE00000340658_SE_7"
    assert se_b["SubexonIDCluster"] == "ENSGGOE00000120742_SE_8"

    assert se_a["TranscriptIDCluster"] == "ENSGGOT00000065893"
    assert se_b["TranscriptIDCluster"] == "ENSGGOT00000065893"

    assert se_a['Strand'] == -1
    assert se_b['Strand'] == -1

    assert se_a["SubexonCodingEnd"] < se_b["SubexonCodingStart"]

    assert se_b["SubexonCodingStart"] == 62034360
    assert se_b["SubexonCodingEnd"] == 62034205

    assert se_a["SubexonCodingStart"] == 62034204
    assert se_a["SubexonCodingEnd"] == 62034178

    assert se_a['StartPhase'] == 1
    assert se_a['EndPhase'] == 1

    assert se_b['StartPhase'] == 1
    assert se_b['EndPhase'] == 1

    assert str(se_b['SubexonSequence']) == (
        "ATATTGATCAGTGGAATAAAGTTATTGAACAGCTTGGAACACCATGTCCTGAATTCATGA"
        "AGAAGCTGCAACCAACAGTAAGGACTTACGTTGAAAACAGACCTAAATATGCTGGATATA"
        "GCTTTGAGAAACTCTTCCCTGATGTCCTTTTCCCAG")
    assert str(se_a['SubexonSequence']) == "CTGACTCAGAACACAACAAACTTAAAG"


def test_subexon_phases_and_coordinates(clustered_trx_data):
    data = clustered_trx_data('MAPK8_all')

    subexon_table = subexons.create_subexon_table(data)

    # ENSABRE00000111150 has two sub-exons: se_a and se_b

    # se_a : AVMSR
    # se_b : SKRDNNFYSVEIGDSTFTVLKRYQNLKPIGSGAQGIVC

    # se_a :
    # \/ 2379656 (start phase: 1)
    # TTGCCGTCATGAGCAG
    #                /\ 2379641 (phase: 2)
    # se_b :
    # \/ 2379640 (phase: 2)
    # AAGCAAGCGCGACAACAATTTCTACAGTGTTGAAATTGGAGATT
    # CTACGTTCACTGTGTTGAAACGGTATCAAAACTTGAAACCAATAGGATCAGGAGCACAGG
    # GGATAGTATG
    #         /\ 2379527 (end phase: 2)

    se_a = subexon_table.loc[1, :]
    se_b = subexon_table.loc[3, :]

    assert se_a["SubexonID"] == "ENSABRE00000111150_SE_14"
    assert se_b["SubexonID"] == "ENSABRE00000111150_SE_13"

    assert se_a["SubexonIDCluster"] == "ENSABRE00000111150_SE_14"
    assert se_b[
        "SubexonIDCluster"] == "ENSABRE00000111150_SE_13/ENSABRE00000111857_SE_13"

    assert se_a["TranscriptIDCluster"] == "ENSABRT00000020442"
    assert se_b["TranscriptIDCluster"] == "ENSABRT00000020442"

    assert se_a['Strand'] == -1
    assert se_b['Strand'] == -1

    assert se_a["SubexonCodingStart"] > se_b["SubexonCodingEnd"]

    assert se_a["SubexonCodingStart"] == 2379656
    assert se_a["SubexonCodingEnd"] == 2379641

    assert se_b["SubexonCodingStart"] == 2379640
    assert se_b["SubexonCodingEnd"] == 2379527

    assert se_a['StartPhase'] == 1
    assert se_a['EndPhase'] == 2

    assert se_b['StartPhase'] == 2
    assert se_b['EndPhase'] == 2

    assert str(se_b['SubexonSequence']) == (
        "AAGCAAGCGCGACAACAATTTCTACAGTGTTGAAATTGGAGATT"
        "CTACGTTCACTGTGTTGAAACGGTATCAAAACTTGAAACCAATAGGATCAGGAGCACAGG"
        "GGATAGTATG")
    assert str(se_a['SubexonSequence']) == "TTGCCGTCATGAGCAG"

    assert str(se_b['SubexonProteinSequence']
               ) == "SKRDNNFYSVEIGDSTFTVLKRYQNLKPIGSGAQGIVC"
    assert str(se_a['SubexonProteinSequence']) == "AVMSR"
