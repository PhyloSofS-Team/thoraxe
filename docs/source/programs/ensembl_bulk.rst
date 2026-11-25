ensembl_bulk
============

`ensembl_bulk` prepares a folder of Ensembl bulk files that `transcript_query`
reads with `--bulk-dir`, instead of querying the Ensembl REST API and BioMart
gene by gene. All the inputs then come from one Ensembl release, whatever the
date of the run, and building the input folder of a gene takes seconds:

::

    ensembl_bulk download ensembl_115 --release 115
    ensembl_bulk index ensembl_115
    transcript_query ENSG00000107643 --bulk-dir ensembl_115 --orthology 1:n

`download` fetches the Compara gene trees and homologies, and, for each
species, the GTF, the unmasked toplevel genome and part of the core database
dump from the Ensembl FTP site, and one exon table from the BioMart archive
site of the release. By default the species are the query species and every
species with an ortholog of it, i.e. those `transcript_query` uses without
`--specieslist`: 199 species and about 120 GB for human in Ensembl 115.
`--specieslist` restricts the download to some species. Every file is checked
against the Ensembl checksums, and running the command again resumes an
interrupted download. `index` then builds the indexes that `transcript_query`
reads. With the decompressed genomes, a species takes about 3 GB (about
0.55 TB for the 199 species).
``scripts/ensembl_bulk.sbatch`` runs both steps as one SLURM job.

`transcript_query --bulk-dir` writes the same files as `transcript_query`,
with two expected differences: the order of sibling subtrees in `tree.nh` can
differ, as it does between two REST queries, and `sequences.fasta` can contain
some exons of overlapping genes that the REST API omits, which `thoraxe` does
not read. Without `--specieslist`, the species of the bulk folder are used.

.. argparse::
   :ref: thoraxe.bulk_ensembl.cli.parse_command_line
   :prog: ensembl_bulk


.. note::
   The bulk mode was validated on Ensembl 115 by running `thoraxe` on 19
   species from REST and from bulk input folders, for the 23,262 human genes of
   a production run. Outputs were equivalent for 19,304 genes. Both sources
   found no 1:n ortholog in the 19 species for 61 genes and no ortholog at all
   for 3,855, and `thoraxe` failed on both for 9. Each of the 33 remaining
   genes had a defect in its REST download (corrupted or missing exon sequences
   of an ortholog, orthologs without exons, empty exon table) or made `thoraxe`
   fail on the extra species of the REST folder.
