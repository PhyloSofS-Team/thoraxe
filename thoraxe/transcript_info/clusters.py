"""
clusters: Helper module for work with sets as clusters.

A cluster is a set of elements and we can have list of clusters.
"""

import logging


def fill_clusters(clusters, element_i, element_j):
    """
    Fill the list of clusters (sets) with two elements from the same cluster.

    >>> clusters = []
    >>> fill_clusters(clusters, 1, 10)
    [{1, 10}]
    >>> fill_clusters(clusters, 10, 100)
    [{1, 10, 100}]
    >>> fill_clusters(clusters, 20, 30)
    [{1, 10, 100}, {20, 30}]
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
    Return the set to delete given a list of sets.

    It keeps the first element of the set.

    >>> set_to_delete([{1, 2, 3}, {4, 5}]) # it keeps 1 and 4
    {2, 3, 5}
    """
    return {element for group in set_list for element in list(group)[1:]}


def inform_about_deletions(todelete, message):
    """
    Warning about elements that are going to be deleted.

    `inform_about_deletions({1, 2, 3}, "Identical sequences were found:")` will inform:
    ::

        WARNING:root:Identical sequences were found:
        WARNING:root:deleting 1
        WARNING:root:deleting 2
        WARNING:root:deleting 3
    """
    if len(todelete) >= 1:
        logging.warning(message)
        for element in todelete:
            logging.warning('deleting %s', element)


def cluster2str(cluster, delim='/', item2str=str):
    """
    Take a cluster (set) and return a string representation.

    >>> cluster2str({1, 10, 100})
    '1/10/100'
    >>> cluster2str({1, 10, 100}, delim=';')
    '1;10;100'
    >>> names = {1 : 'a', 10 : 'b'}
    >>> cluster2str({1, 10, 100}, item2str=lambda x: names.get(x, str(x)))
    'a/b/100'
    """
    return delim.join(item2str(element) for element in sorted(cluster))
