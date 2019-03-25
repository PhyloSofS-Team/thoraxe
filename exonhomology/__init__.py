"""
exonhomology: Pipeline to disentangle homology relationships between exons.

It uses the other submodules to offer a higher level interface.
"""

from exonhomology import utils
from exonhomology import transcript_query
from exonhomology import transcript_info
from exonhomology import subexons

from exonhomology.exonhomology import *
