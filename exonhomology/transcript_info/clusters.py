"""
clusters: Helper module for work with sets as clusters.
=======================================================

A cluster is a set of elements and we can have list of clusters.
"""


def fill_clusters(clusters, element_i, element_j):
    """
    Fill the list of clusters (sets) with two elements from the same cluster.

    >>> clusters = []
    >>> fill_clusters(clusters, 1, 10)
    [set([1, 10])]
    >>> fill_clusters(clusters, 10, 100)
    [set([1, 10, 100])]
    >>> fill_clusters(clusters, 20, 30)
    [set([1, 10, 100]), set([20, 30])]
    """
    if not clusters:  # i.e. len(clusters) == 0
        clusters.append({element_i, element_j})
    else:
        found = False
        for cluster in clusters:
            if element_i in cluster:
                cluster.add(element_j)
                found = True
            elif element_j in cluster:
                cluster.add(element_i)
                found = True
        if not found:
            clusters.append({element_i, element_j})
    return clusters


def set_to_delete(set_list):
    """
    It returns the set to delete given a list of sets.
    It keeps the first element of the set.

    >>> set_to_delete([set([1, 2, 3]), set([4, 5])]) # it keeps 1 and 4
    set([2, 3, 5])
    """
    return {element for group in set_list for element in list(group)[1:]}
