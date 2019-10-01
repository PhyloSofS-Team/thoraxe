Output description
==================


*ThorAxe* has two command-line programs, `transcript_query` and `thoraxe`. The
first downloads from Ensembl_ the inputs for the second.


transcript_query
----------------

::
    transcript_query PAX6 -l homo_sapiens,mus_musculus,xenopus_tropicalis,danio_rerio

Running the previous line is going to create a `PAX6` folder with an `Ensembl`
subfolder. That last folder contains four files with information downloaded
from *Ensembl*. In this case, the information corresponds to the transcripts
and exons from the human *PAX6* gene and all their one-to-one orthologs in the
indicated species:

::

    PAX6
    └── Ensembl
        ├── exonstable.tsv
        ├── sequences.fasta
        ├── tree.nh
        └── tsl.csv


thoraxe
-------

::

    PAX6
    ├── Ensembl
    │   ├── exonstable.tsv
    │   ├── sequences.fasta
    │   ├── tree.nh
    │   └── tsl.csv
    └── thoraxe
        ├── _intermediate
        │   ├── chimeric_alignment_1.fasta
        │   ⋮
        │   ├── cluster_data.js
        │   ├── cluster_plots.html
        │   ├── gene_ids_1.txt
        │   ⋮
        │   ├── msa_matrix_1.txt
        │   ⋮
        │   ├── subexon_table.csv
        │   ├── subexon_table_1.csv
        │   ⋮
        │   └── transcript_table.csv
        ├── homologous_exon_table.csv
        ├── msa
        │   ├── msa_homologous_exon_10_0.fasta
        │   ⋮
        │   └── msa_homologous_exon_9_0.fasta
        ├── phylosofs
        │   ├── homologous_exons.tsv
        │   ├── transcripts.pir
        │   ├── transcripts.txt
        │   └── tree.nh
        └── splice_graph.gml


For each protein, the `Ensembl` folder has the data downloaded from Ensembl
with `transcript_query`. The `thoraxe` folder as the raw output from `thoraxe`.


thoraxe
~~~~~~~

The `homologous_exon_table.csv` file has the tidy (denormalized) output data for
homologous exons. Each row is an observation of a homologous exon in a
particular transcript. Rows are in order, following the homologous exon rank.
That means that the concatenation of the homologous exon sequence for each
transcript gives the protein isoform sequence.

The `splice_graph.gml` is the splice graph of homologous exons in all the
transcripts/species using the human-readable Graph Modelling Language (GML)
format. Each node and edge has conservation information, where conservation
means the fraction of species showing that particular node (homologous exon) or
connection.

The `msa` folder has a multiple sequence alignment for each homologous exon, e.g.
`msa_homologous_exon_10_0.fasta`. This MSAs can be easily used as a seed to
look for homologous sequences in different databases by using
[hmmsearch](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch).


phylosofs
~~~~~~~~~

The `phylosofs` folder has the needed inputs for the structural and mollecular
modelling pipelines of *PhyloSofS*. In particular, these files use a single
unicode character to represent each homologous exon. The mapping between the id
for *PhyloSofS* and the one of *ThorAxe* (*ExonCluster*_*ChimericBlock*) is in
the `homologous_exons.tsv` file.

The `transcripts.pir` has the annotated sequence using the PIR format.

The list of transcripts for each gene is in `transcripts.txt`


Intermediate outputs
~~~~~~~~~~~~~~~~~~~~

The `_intermediate` folder has intermediate files and states from the *ThorAxe*
pipeline. In particular, `cluster_plots.html` has the interactive plots of the
chimeric alignments, with the possibility to show the `constitutive` value
calculated for each subexon as described in previous issues of this repo.



.. _Ensembl: https://www.ensembl.org/index.html
