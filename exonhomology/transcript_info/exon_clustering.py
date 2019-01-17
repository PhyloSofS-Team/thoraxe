from skbio import Protein
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment._pairwise import blosum50


def exon_clustering(trx_data):
    """
    Cluster exons based on their sequence identity after local alignment.
    """
    local_pairwise_align_ssw(
        Protein("SS"),
        Protein("SS"),
        gap_open_penalty=10,
        gap_extend_penalty=1,
        substitution_matrix=blosum50)
