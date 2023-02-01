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
        ├── ensembl_version.csv
        ├── exonstable.tsv
        ├── sequences.fasta
        ├── tree.nh
        └── tsl.csv

The `ensembl_version.csv` file contains a table with information about the
download date and hour together with the *Ensembl Genomes* version.

`sequences.fasta` contains the CDS sequence of each downloaded exon. You can
find the used genome assembly in the sequence name/description. For example:

::

    >homo_sapiens:MAPK8 ENSE00003837028 chromosome:GRCh38:10:48306753:48306821:1
    CGGCGACCACCCCGGACGGCCCCTGTCCCCGCTGGCGGGCTTCCCTGTCGCCGTTCGCTGCGCTGCCGG

It's the entry for **MAPK8** exon `ENSE00003837028` from chromosome **10** in
the assembly **GRCh38**.

thoraxe
-------

We can run `thoraxe` in the `PAX6` folder that `transcript_query` created with
the *Ensembl* data:

::

    thoraxe -i PAX6 --phylosofs --plot_chimerics

*ThorAxe* is going to look for the `Ensembl` folder with the input data and is
going to create a `thoraxe` folder containing the final and intermediate
results of the study. The `--phylosofs` and `--plot_chimerics` arguments are
optionals.

::

    PAX6
    ├── Ensembl
    │   ├── ensembl_version.csv
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
        ├── ases_table.csv
        ├── msa
        │   ├── msa_s_exon_10_0.fasta
        │   ⋮
        │   └── msa_s_exon_9_0.fasta
        ├── path_table.csv
        ├── phylosofs
        │   ├── s_exons.tsv
        │   ├── transcripts.pir
        │   ├── transcripts.txt
        │   └── tree.nh
        ├── s_exon_table.csv
        └── splice_graph.gml


The `msa` folder has a multiple sequence alignment for each s-exon, e.g.
`msa_s_exon_10_0.fasta`. This MSAs can be easily used as a seed to
look for homologous sequences in different databases by using hmmsearch_.


Tables
~~~~~~

*ThorAxe* outputs three tidy (denormalized) tables rather than a normalized
database: `s_exon_table.csv`, `path_table.csv` and `ases_table.csv`.  These
tables are comma-separated files easily read by R_, Python_ (using pandas_) or
Julia_ (using CSV_ and DataFrames_) among other software tools. When a cell
should include a list of elements, *ThorAxe* uses by default the character `/`
as the delimiter between the list elements. The stored list can easily
re-create by using `strsplit` in *R*, the `split` string method in *Python* or
the `split` function in *Julia*.


S-exon table
''''''''''''

The `s_exon_table.csv` file has the tidy (denormalized) output data for
the **orthologous exonic regions (s-exons)**. Each row is an observation of an
s-exon in a particular transcript. Rows are in order, following the s-exon rank.
That means that the concatenation of the s-exon sequences for each
transcript gives the protein isoform sequence.


Path table
''''''''''

The `path_table.csv` contains information about each transcript as a path in
the splice graph. The first row indicates the *canonical path*. The canonical
path is used to detect alternative splicing events automatically, and it is
the one present in the largest number of genes/species with the highest
minimal value of `transcript_weighted_conservation`, and it is the longest. The
alternative splicing events are stored in the `ases_table.csv` file where the
canonical path and the alternative are indicated, together with conservation
(number of genes showing the path) information.


Table of Alternative Splicing Events
''''''''''''''''''''''''''''''''''''

ThorAxe automatically detects alternative splicing events by comparing each
transcript/path to the canonical path in the splice graph. The `ases_table.csv`
file stores these detected events. The table identifies each event using the
first two columns: `CanonicalPath` and `AlternativePath`. These columns store
the list of nodes that forms the alternative subpaths in the canonical and
*alternative* (non-canonical) path. First, *ThorAxe* classifies each event,
respect to the canonical path, as `insertion`, `deletion` ,
`fully_alternative`, `alternative_start`, `alternative_end` or simply
`alternative`.  The `ASE` column stores this classification. The events,
excepting insertions and deletions, are further classified according to their
number of mutually exclusive s-exons pairs into: `mutually_exclusive`,
`partially_mutually_exclusive` and not mutually exclusive (empty cell). The
`MutualExclusivity` column stores this classification. The
`MutualExclusiveCanonical` and `MutualExclusiveAlternative` columns stores two
paired list indicating the mutually exclusive s-exon pairs. For example, an
event having `15_0/15_1` in `MutualExclusiveCanonical` and `7_5/7_5` in
`MutualExclusiveAlternative`, means that s-exon `15_0` is mutually exclusive
with s-exon `7_5` and that `15_1` is also mutually exclusive with `7_5`. Here,
mutually exclusive means that there is no transcript in the path table showing
both s-exons.


Splice graph
~~~~~~~~~~~~

*ThorAxe* splice graph as *s-exons*, rather than exon, as nodes. Using *s-exon*
allows us to represent transcripts from different species in the same splice
graph. As a consequence, the edges indicate the connectivity between these
exonic regions and not necessarily junctions between genomic exons. *ThorAxe*
splice graph has conservation information, where conservation means the
fraction or transcript weighted fraction of species showing a particular node
(s-exon) or connection. The `splice_graph.gml` file stores this splice graph
of s-exons using the human-readable *Graph Modelling Language* (GML) format.
*GML* is a rich format, and we used it to store useful metadata for nodes and
edges:


================================ ===============================================
            Metadata                              Description
================================ ===============================================
label                            The s-exon ID (only for nodes).
transcript_fraction              It is the number of transcripts showing that
                                 node/edge over the total number of transcripts,
                                 taking into account all the genes.
conservation                     Fraction of genes showing that node/edge.
transcript_weighted_conservation As `conservation` but each gene is weighted
                                 using the fraction of transcripts in that
                                 gene showing that edge (only for edges).
genes                            List of genes, separated by commas,
                                 showing that node/edge.
transcripts                      List of transcripts, separated by commas,
                                 showing that node/edge.
================================ ===============================================


.. note::
    If *ThorAxe's* `transcript_query` is run with `--orthology 1:1` (default),
    the number of genes is identical to the number of species in the dataset.


phylosofs
~~~~~~~~~

The `phylosofs` folder has the needed inputs for the structural and molecular
modelling pipelines of PhyloSofS_. It is only generated when the `--phylosofs`
optional argument is used.
These PhyloSofS's files use a single unicode character to represent each s-exon.
The mapping between the id for *PhyloSofS* and the one of *ThorAxe*
(*ExonCluster_ChimericBlock*) is in the `s_exons.tsv` file.

The `transcripts.pir` has the annotated sequence using the PIR format.

The list of transcripts for each gene is in `transcripts.txt`

Intermediate outputs
~~~~~~~~~~~~~~~~~~~~

The `_intermediate` folder has intermediate files and states from the *ThorAxe*
pipeline. In particular, `cluster_plots.html` has the interactive plots of the
chimeric alignments, with the possibility to show the `constitutive` value
calculated for each subexon. This plot is only generated when the
`--plot_chimerics` optional argument is indicated.

.. warning::
    `cluster_plots.html` loads the `cluster_data.js` file that includes
    all the information needed to plot all the generated multiple sequence alignments
    using *Plotly*. It can freeze the tab when a lot of sequences/species are used or
    if the computer hasn't enough resources available.


.. _Ensembl: https://www.ensembl.org/index.html
.. _hmmsearch: //www.ebi.ac.uk/Tools/hmmer/search/hmmsearch
.. _PhyloSofS: https://github.com/PhyloSofS-Team/PhyloSofS
.. _R: https://www.r-project.org/
.. _Python: https://www.python.org/
.. _Julia: https://julialang.org/
.. _pandas: https://pandas.pydata.org/
.. _CSV: https://github.com/JuliaData/CSV.jl
.. _DataFrames: https://github.com/JuliaData/DataFrames.jl
