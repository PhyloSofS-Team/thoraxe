"""
bulk_ensembl: build transcript_query's input folders from Ensembl bulk files.

`transcript_query` asks the Ensembl REST API and BioMart, gene by gene, for the
files that `thoraxe` reads from ``<gene>/Ensembl/``. This subpackage builds the
same files from the release-wide files that Ensembl publishes on its FTP site
and from one BioMart dump per species, downloaded once:

==================== ======================================= ====================
file                 online source (transcript_query)         bulk source
==================== ======================================= ====================
tree.nh              REST /genetree/member/id                 Compara EMF dump
orthologs            REST /homology/id                        Compara homology TSV
tsl.csv              REST /overlap/id?feature=transcript      core MySQL dump
exonstable.tsv       BioMart, one query per gene              BioMart, one query
                                                              per species
sequences.fasta      REST /overlap/id?feature=exon and        genome FASTA + GTF
                     /sequence/id
ensembl_version.csv  REST /info/eg_version                    release number
==================== ======================================= ====================

The bulk folder is prepared once with ``ensembl_bulk download`` and
``ensembl_bulk index``, then ``transcript_query --bulk-dir <folder>`` reads it.
Every query is then local, so the inputs of a run do not depend on the state of
the Ensembl servers at download time, and they always come from one release.
"""

from thoraxe.bulk_ensembl.layout import BulkDir
