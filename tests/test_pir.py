import os
import pytest
from Bio import SeqIO

@pytest.fixture(scope='module')
def mapk8(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'data')
    mapk8_dir = os.path.join(data_dir, 'MAPK8_output', 'thoraxe')
    return {
        'transcripts': os.path.join(mapk8_dir, 'phylosofs', 'transcripts.pir')
    }


@pytest.fixture(scope='module')
def polr3b(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'data')
    polr3b_dir = os.path.join(data_dir, 'POLR3B_output', 'thoraxe')
    return {
        'transcripts': os.path.join(polr3b_dir, 'phylosofs', 'transcripts.pir')
    }


def test_shared_codon(polr3b):
    """
    Three sub-exons sharing a codon in reverse strand:

    ```
    ENSSSCE00000329694_SE_19: ATGCATTGTGGGTTACTGGAGGAACCCGACATGGATTCCACAG
    ENSSSCE00000376311_SE_18: A
    ENSSSCE00000376311_SE_17: GAGCTGGATTGAGAGATGTCTCAACGAAAGTGAAAACAAACGTTATTCCAGCCA...
    ```
    """
    seqs = SeqIO.to_dict(SeqIO.parse(polr3b['transcripts'], "pir"))
    keys = [key for key in seqs.keys() if "ENSSSCT00000000177" in key]
    assert keys
    assert seqs[keys[0]].seq == """
    MHCGLLEEPDMDSTESWIERCLNESENKRYSSHTSLGNVSNDENEEKENNRASKPHSTPA
    TLQWLEENYEIAEGVCIPRSALYMHYLDFCEKNDTQPVNAASFGKIIRQQFPQLTTRRLG
    TRGQSKYHYYGIAVKESSQYYDVMYSKKGAAWVSETGKKEVSKQTVAYSPRSKLGTLLPE
    FPNVKDLNLPASLPEEKVSTFIMMYRTHCQRILDTVIRANFDEVQSFLLHFWQGMPPHML
    PVLGSSTVVNIVGVCDSILYKAISGVLMPTVLQALPDSLTQVIRKFAKQLDEWLKVALHD
    LPENLRNIKFELSRRFSQILRRQTSLNHLCQASRTVIHSADITFQMLEDWRNVDLNSITK
    QTLYTMEDTRDEHRKLIIQLYQEFDHLLEEQSPIESYIEWLDTMVDRCVVKVAAKRQGSL
    KKVAQQFLLMWSCFGTRVIRDMTLHSAPSFGSFHLIHLMFDDYVLYLLESLHCQERASEL
    MRAMKGEGSTADVREDIILAEAAPPTPSPVPSFSPAKSATSVEVPPPSSPVSNPSPEYTG
    LSTTGAMQSYTWSLTYTVTTAAGSPAENSQQLPCMRSTHVPSSSVTHRIPVYPHREEHGY
    TGSYNYGSYGNQHPPPMQSQYPALPHDTAISGPLHYSPYHRSSAQYPFNSPTSRMEPCLM
    SSTPRLHPTPVTPRWPEVPTANTCYTSPSVHSTRYGNSSDMYTPLTTRRNSEYEHMQHFP
    GFAYINGEASTGWAK
    """.replace("\n", "").replace(" ", "")
    

def test_polr3b_human(polr3b):
    # Human : forward strand and shared codon
    seqs = SeqIO.to_dict(SeqIO.parse(polr3b['transcripts'], "pir"))
    keys = [key for key in seqs.keys() if "ENST00000228347" in key]
    assert keys
    assert seqs[keys[0]].seq == """
    MDVLAEEFGNLTPEQLAAPIPTVEEKWRLLPAFLKVKGLVKQHIDSFNYFINVEIKKIMK
    ANEKVTSDADPMWYLKYLNIYVGLPDVEESFNVTRPVSPHECRLRDMTYSAPITVDIEYT
    RGSQRIIRNALPIGRMPIMLRSSNCVLTGKTPAEFAKLNECPLDPGGYFIVKGVEKVILI
    QEQLSKNRIIVEADRKGAVGASVTSSTHEKKSRTNMAVKQGRFYLRHNTLSEDIPIVIIF
    KAMGVESDQEIVQMIGTEEHVMAAFGPSLEECQKAQIFTQMQALKYIGNKVRRQRMWGGG
    PKKTKIEEARELLASTILTHVPVKEFNFRAKCIYTAVMVRRVILAQGDNKVDDRDYYGNK
    RLELAGQLLSLLFEDLFKKFNSEMKKIADQVIPKQRAAQFDVVKHMRQDQITNGMVNAIS
    TGNWSLKRFKMDRQGVTQVLSRLSYISALGMMTRISSQFEKTRKVSGPRSLQPSQWGMLC
    PSDTPEGEACGLVKNLALMTHITTDMEDGPIVKLASNLGVEDVNLLCGEELSYPNVFLVF
    LNGNILGVIRDHKKLVNTFRLMRRAGYINEFVSISTNLTDRCVYISSDGGRLCRPYIIVK
    KQKPAVTNKHMEELAQGYRNFEDFLHESLVEYLDVNEENDCNIALYEHTINKDTTHLEIE
    PFTLLGVCAGLIPYPHHNQSPRNTYQCAMGKQAMGTIGYNQRNRIDTLMYLLAYPQKPMV
    KTKTIELIEFEKLPAGQNATVAVMSYSGYDIEDALVLNKASLDRGFGRCLVYKNAKCTLK
    RYTNQTFDKVMGPMLDAATRKPIWRHEILDADGICSPGEKVENKQVLVNKSMPTVTQIPL
    EGSNVPQQPQYKDVPITYKGATDSYIEKVMISSNAEDAFLIKMLLRQTRRPEIGDKFSSR
    HGQKGVCGLIVPQEDMPFCDSGICPDIIMNPHGFPSRMTVGKLIELLAGKAGVLDGRFHY
    GTAFGGSKVKDVCEDLVRHGYNYLGKDYVTSGITGEPLEAYIYFGPVYYQKLKHMVLDKM
    HARARGPRAVLTRQPTEGRSRDGGLRLGEMERDCLIGYGASMLLLERLMISSDAFEVDVC
    GQCGLLGYSGWCHYCKSSCHVSSLRIPYACKLLFQELQSMNIIPRLKLSKYNE
    """.replace("\n", "").replace(" ", "")


def test_polr3b_consistency(polr3b):
    seqs = SeqIO.to_dict(SeqIO.parse(polr3b['transcripts'], "pir"))
    for (seq_id, seq) in seqs.items():
        assert len(seq.seq) == len(seq.description)

def test_polr3b_consistency(mapk8):
    seqs = SeqIO.to_dict(SeqIO.parse(mapk8['transcripts'], "pir"))
    for (seq_id, seq) in seqs.items():
        assert len(seq.seq) == len(seq.description)
