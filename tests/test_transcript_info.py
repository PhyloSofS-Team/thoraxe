import os
import pytest
import numpy as np
from thoraxe import transcript_info


@pytest.fixture(scope='module')
def datadir(request):
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    return os.path.join(test_dir, 'data')


@pytest.fixture(scope='module')
def mapk8(datadir):
    mapk8_dir = os.path.join(datadir, 'MAPK8', 'Ensembl')
    return {
        'tsl': os.path.join(mapk8_dir, 'tsl.csv'),
        'exontable': os.path.join(mapk8_dir, 'exonstable.tsv'),
        'seqs': os.path.join(mapk8_dir, 'sequences.fasta')
    }


@pytest.fixture(scope='module')
def camk2a(datadir):
    camk2a_dir = os.path.join(datadir, 'CAMK2A', 'Ensembl')
    return {
        'tsl': os.path.join(camk2a_dir, 'tsl.csv'),
        'exontable': os.path.join(camk2a_dir, 'exonstable.tsv'),
        'seqs': os.path.join(camk2a_dir, 'sequences.fasta')
    }


@pytest.fixture(scope='module')
def grin1(datadir):
    grin1_dir = os.path.join(datadir, 'GRIN1', 'Ensembl')
    return {
        'tsl': os.path.join(grin1_dir, 'tsl.csv'),
        'exontable': os.path.join(grin1_dir, 'exonstable.tsv'),
        'seqs': os.path.join(grin1_dir, 'sequences.fasta')
    }


def test_read_transcript_info(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'], mapk8['exontable'], mapk8['seqs'])

    assert trx_data.loc[trx_data['TranscriptStableID'] ==
                        'ENST00000374179', 'Flags'].unique(
                        )[0] == '1 (assigned to previous version 7)'

    # Only h. sapiens & m. musculus have TSL information:
    assert all(trx_data.Species.unique() == ['homo_sapiens', 'mus_musculus'])


def test_remove_na(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'], mapk8['exontable'], mapk8['seqs'], remove_na=False)

    # I keep other species, not only h. sapiens & m. musculus:
    assert len(trx_data.Species.unique()) > 2

    # remove_na doesn't interfere with selecting the correct biotype
    assert sum(value not in ['Protein coding', 'protein_coding']
               for value in trx_data.Biotype.unique()) == 0

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert sum(trx_data['TranscriptStableID'] == 'ENSRNOT00000083933') == 0


def test_species_list(mapk8):
    species_list = ['bos_taurus', 'danio_rerio', 'homo_sapiens']

    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs'],
        remove_na=False,
        species_list=species_list)

    assert sorted(trx_data.Species.unique()) == species_list


def test_keep_badquality_sequences(mapk8):
    trx_data = transcript_info.read_transcript_info(
        mapk8['tsl'],
        mapk8['exontable'],
        mapk8['seqs'],
        remove_na=False,
        remove_badquality=False)

    # ENSRNOT00000083933 has Xs in its sequence: ...VILGMGYKENGQXVXHVQRGLICC*
    assert sum(trx_data['TranscriptStableID'] == 'ENSRNOT00000083933') == 5


def test_non_coding_exons_camk2a(camk2a):
    trx_data = transcript_info.read_transcript_info(
        camk2a['tsl'], camk2a['exontable'], camk2a['seqs'], remove_na=False)

    # The two first exons of ENSSSCT00000052397 are non-coding
    assert ''.join(
        str(exon)
        for exon in trx_data.loc[trx_data['TranscriptStableID'] ==
                                 'ENSSSCT00000052397', 'ExonProteinSequence']
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
        str(exon) for exon in
        trx_data.loc[trx_data['TranscriptStableIDCluster'].
                     map(lambda ids: 'ENSMUST00000099506' in ids.split('/')
                         ), 'ExonProteinSequence']
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
    assert all(clustered['ExonStableID'] == trx_data['ExonStableID'])

    # Exon can not have more than one cluster
    assert all(
        clustered.groupby(
            'ExonStableID').apply(lambda df: len(df['Cluster'].unique()) == 1))

    # Non-clustered exons have Cluster == 0 and QueryExon == ''
    assert all(row.QueryExon == '' for row in clustered.itertuples()
               if row.Cluster == 0)
    assert all(row.Cluster == 0 for row in clustered.itertuples()
               if row.QueryExon == '')

    # Sequences with less than 4 residues are non-clustered by default
    assert all(
        clustered.loc[clustered['ExonProteinSequence'].map(len) < 4, 'Cluster']
        == 0)

    for _, group in clustered.groupby('Cluster'):
        nans = np.isnan(group['PercentIdentity'])
        # There is a nan in PercentIdentity when a sequence initialize its own
        # cluster, so the QueryExon and the ExonStableID should be the same
        assert np.all(
            group.loc[nans, 'QueryExon'] == group.loc[nans, 'ExonStableID'])
        # Also, if there are more exons, the exons with nan should be the
        # QueryExon of other exon in the cluster. It can not be alone.
        if len(group) > 1:
            assert np.all([
                exon in group.loc[np.logical_not(nans), 'QueryExon'].values
                for exon in group.loc[nans, 'QueryExon'].unique()
            ])

            # The aligned seq should be in the exon
            subset = group[np.logical_not(nans)]
            for _, exon in subset.groupby('ExonStableIDCluster'):
                assert np.any([
                    row.AlignedTarget.replace('-',
                                              '') in row.ExonProteinSequence
                    for row in exon.itertuples()
                ])
