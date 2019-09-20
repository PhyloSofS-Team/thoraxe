"""
species: Util functions to get the species list.
"""

import os
import re


def check_species_name(species_name):
    """
    Return True if the species_name has the correct format.

    The expected name is binomial but it can also include the subspecies name
    (trinomial). If the species_name doesn't conform the expected format an
    error is raised.

    >>> check_species_name('homo_sapiens')
    True
    >>> check_species_name('colobus_angolensis_palliatus')
    True
    >>> check_species_name('cricetulus_griseus_chok1gshd')
    True
    """
    result = re.match('^[a-z]+_[a-z]+(_[0-9a-z]+)?$', species_name)
    if result is None:
        raise ValueError(
            'The species name should be the binomial/trinomial name lowercased'
            ' with spaces replaced by underscores e.g. Homo sapiens should be '
            'homo_sapiens. {} do not conform the format.'.format(species_name))
    return True


def _read_species_list(file):
    """Return the species list from the file."""
    species = []
    with open(file, 'r') as stream:
        for line in stream:
            species_name = line.strip()
            if check_species_name(species_name):
                species.append(species_name)
    return species


def get_species_list(specieslist):
    """
    Return a list of species names (str).

    >>> get_species_list('')
    >>> get_species_list('homo_sapiens,mus_musculus')
    ['homo_sapiens', 'mus_musculus']
    """
    if specieslist == '':
        return None

    fields = specieslist.split(',')

    if len(fields) == 1:
        filename = fields[0]
        if os.path.isfile(filename):
            return _read_species_list(filename)
        raise ValueError('{} is not a file.'.format(filename))

    for species_name in fields:
        check_species_name(species_name)

    return fields
