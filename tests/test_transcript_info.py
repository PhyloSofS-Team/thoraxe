import os
import pytest
from exonhomology import transcript_info


@pytest.fixture(scope='module')
def datadir(request):
    return os.path.join('.', 'tests', 'data')


@pytest.fixture(scope='module')
def mapk8(datadir):
    mapk8_dir = os.path.join(datadir, 'MAPK8_ENSG00000107643')
    return {
        'tsl':
        os.path.join(mapk8_dir, 'TSL', 'MAPK8_ENSG00000107643_TSL.csv'),
        'exontable':
        os.path.join(mapk8_dir, 'TablesExons',
                     'MAPK8_ENSG00000107643_exonstable.tsv'),
        'seqs':
        os.path.join(mapk8_dir, 'Sequences', 'MAPK8_ENSG00000107643.fasta')
    }


@pytest.fixture(scope='module')
def camk2a(datadir):
    camk2a_dir = os.path.join(datadir, 'CAMK2A_ENSG00000070808')
    return {
        'tsl':
        os.path.join(camk2a_dir, 'TSL', 'CAMK2A_ENSG00000070808_TSL.csv'),
        'exontable':
        os.path.join(camk2a_dir, 'TablesExons',
                     'CAMK2A_ENSG00000070808_exonstable.tsv'),
        'seqs':
        os.path.join(camk2a_dir, 'Sequences', 'CAMK2A_ENSG00000070808.fasta')
    }


@pytest.fixture(scope='module')
def grin1(datadir):
    grin1_dir = os.path.join(datadir, 'GRIN1_ENSG00000169258')
    return {
        'tsl':
        os.path.join(grin1_dir, 'TSL', 'GRIN1_ENSG00000169258_TSL.csv'),
        'exontable':
        os.path.join(grin1_dir, 'TablesExons',
                     'GRIN1_ENSG00000169258_exonstable.tsv'),
        'seqs':
        os.path.join(grin1_dir, 'Sequences', 'GRIN1_ENSG00000169258.fasta')
    }


def test_read_transcript_info(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'], mapk8['exontable'], mapk8['seqs'])

    assert trx_data.loc[trx_data['Transcript stable ID'] == 'ENST00000374179',
                        'Flags'].unique()[
                            0] == '1 (assigned to previous version 7)'

    # Only h. sapiens & m. musculus have TSL information:
    assert all(trx_data.Species.unique() == ['homo_sapiens', 'mus_musculus'])


def test_remove_na(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'], mapk8['exontable'], mapk8['seqs'], remove_na=False)

    # I keep other species, not only h. sapiens & m. musculus:
    assert len(trx_data.Species.unique()) > 2

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert 0 == sum(trx_data['Transcript stable ID'] == 'ENSRNOT00000083933')


def test_keep_badquality_sequences(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs'],
        remove_na=False,
        remove_badquality=False)

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert 5 == sum(trx_data['Transcript stable ID'] == 'ENSRNOT00000083933')


def test_non_coding_exons_camk2a(camk2a):
    trx_data = transcript_info.read_transcript_info(
        camk2a['tsl'], camk2a['exontable'], camk2a['seqs'], remove_na=False)

    # The two first exons of ENSSSCT00000052397 are non-coding
    assert ''.join(
        str(exon)
        for exon in trx_data.loc[trx_data['Transcript stable ID'] ==
                                 'ENSSSCT00000052397', 'Exon protein sequence']
    ) == ('MLLFLALWALVPCLVLLSLYFYSSAGGKSGGNKKNDGVKKRKSSSSVQLMESSESTNTTI'
          'EDEDTKVRKQEIIKVTEQLIEAISNGDFESYTKMCDPGMTAFEPEALGNLVEGLDFHRFY'
          'FENLWSRNSKPVHTTILNPHIHLMGDESACIAYIRITQYLDAGGIPRTAQSEETRVWHRR'
          'DGKWQIVHFHRSGAPSVLPH*')


def test_non_coding_exons_grin1(grin1):
    trx_data = transcript_info.read_transcript_info(
        grin1['tsl'], grin1['exontable'], grin1['seqs'], remove_na=False)

    # The last two exons of ENSMUST00000099506 are non-coding.
    # The unique coding exon has UTR at both ends.
    assert ''.join(
        str(exon)
        for exon in trx_data.loc[trx_data['Transcript stable ID'] ==
                                 'ENSMUST00000099506', 'Exon protein sequence']
    ) == ('MRDCCSSPKAIPAPPRHALDQSLGMDPRHTSSSGAAEGASCSERPAGSLACPSPNCSPLP'
          'ETPRAHGALTSDNSGTTLFGKPEPMSSAEATPTASEIRNPVFSGKMDGNSLKQADSTSTR'
          'KEEAGSLRNEESMLKGKAEPMIYGKGEPGTVGRVDCTASGAENSGSLGKVDMPCSSKVDI'
          'VSPGGDNAGSLRKVETISSGKMDPKTENVMHSRRERPGSTGEGDLVSLRENDMKPPDNTD'
          'SASTKKTDPEFSGKLTPGSSGKTELVSSVTVAPVTSENVNPVCSGGAGPAAVGNSETLSS'
          'VKKDPQLLGKKEAVSSGEGGSVSVRMAETVSARQPEGMFPAKTDSTSSNSTGPSGRADPV'
          'SLRNSELVSPVKPERLSSGQAERVSLVKTETLSSGKEDPRSSRRVDHTTVTGNMQTSQKG'
          'NPESSGKTDLGSSSSGDTRSLGTWGSLSAAKAEVTEGKGDPQPWKKASLPASEKTDPLAS'
          'SKAGSASQGKAETVSPGEVDAMTLGKTVPTSSGKTALVSPGKVDLMTSERAEGIPELQAS'
          'EKGNPVNSTRVDTGATGSTEPKSGVKVITQIPGATSPGKVETPSLQKEQPQLSEKTDPSR'
          'KVDPPTTVEPVSLGKADSASPSPRKAESQTSAKTVPQAPDKATSSLRQSDGTPYSSAQPQ'
          'RDTRSIGSLPEREPSASTSQKDLAAAAAQKSPSAEAAAPPPGPRTRDNFTKAPSWDAGAP'
          'PPREDAGTQAGAQACVSVAVSPMSPQDGAGGPAFSFQAAPRAPSPAPRPPSRRDAGLQVS'
          'LGAAETRSVATGPMTPQAAAPPAVPPVFPEVRVRPGSVLAAALAPQEATEPVRDVSWDEK'
          'GMTWEVYGASMEVEVLGMAIQKHLERQIEEHGRQGAPAPAPPPAVRAGPGRAGSVRTAPA'
          'EGAAKRPPGLFRALLQSVRRPRCCSRAGPTAE*')


def test_exon_clustering(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'], mapk8['exontable'], mapk8['seqs'])

    clustered = transcript_info.exon_clustering(trx_data)

    # Input order doesn't change
    assert all(clustered['Exon stable ID'] == trx_data['Exon stable ID'])

    # Exon can not have more than one cluster
    assert all(
        clustered.groupby('Exon stable ID').apply(
            lambda df: len(df['Cluster'].unique()) == 1))

    # Non-clustered exons have Cluster == 0 and QueryExon == ''
    assert all(row.QueryExon == '' for row in clustered.itertuples()
               if row.Cluster == 0)
    assert all(row.Cluster == 0 for row in clustered.itertuples()
               if row.QueryExon == '')

    # Sequences with less than 4 residues are non-clustered by default
    assert all(clustered.loc[clustered['Exon protein sequence'].map(len) < 4,
                             'Cluster'] == 0)
