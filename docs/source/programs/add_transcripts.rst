add_transcripts
===============

This script should run after `transcript_query` and before `thoraxe`. 
For example, let's see you want to add a *MAPK8* transcript described in 
`user_transcript.csv` to the data download from *Ensembl* for that gene to run 
`thoraxe`; you can do something like:

::

    transcript_query MAPK8
    add_transcripts user_transcript.csv MAPK8/Ensembl
    thoraxe -i MAPK8


.. argparse::
    :ref: thoraxe.add_transcripts.add_transcripts.parse_command_line
    :prog: transcript_query

.. note::
    While this program is called `add_transcripts` you can use it to add genes 
    by adding their transcripts. Use a single input table to add multiple 
    genes and transcripts or run this script multiple times to add a different 
    transcript each time.

.. warning::
    You can not use this program to add single exons unless they contain the 
    complete *CDS*. Otherwise, `thoraxe` would delete the exon and its 
    incomplete transcript.

Input preparation
~~~~~~~~~~~~~~~~~

The input table should be a CSV file with the following columns:

==================== ===========================================================
    Column Name                         Description
==================== ===========================================================
Species              It should be the binomial species name in lowercase and 
                     using underscore instead of space.
GeneID               It should be the Ensembl gene ID, e.g. ENSG00000107643, 
                     rather than the gene name, e.g. MAPK8.
TranscriptID         A string to identify the transcript.
Strand               The Strand should be 1 for a gene in the forward strand 
                     and -1 for one in the reverse strand.
ExonID               A string to identify the exon.
ExonRank             The Exon Rank should be consecutive integer numbers 
                     indicating the order of the exons in the transcript.
ExonRegionStart      ExonRegionStart should be the genomic coordinate of the 
                     first, last if the gene is in the reverse strand, 
                     nucleotide of the NucleotideSequence of the exon. Note 
                     that ExonRegionStart should always be less than 
                     ExonRegionEnd.
ExonRegionEnd        ExonRegionEnd should be the genomic coordinate of the 
                     last, first if the gene is in the reverse strand, 
                     nucleotide of the NucleotideSequence of the exon.
GenomicCodingStart   If the gene is in the forward strand, GenomicCodingStart 
                     should be the first coding nucleotide's genomic 
                     coordinate. Otherwise, it should be the genomic 
                     coordinate of the last coding nucleotide. Note that 
                     GenomicCodingStart should always be less than 
                     GenomicCodingEnd.
GenomicCodingEnd     If the gene is in the forward strand, GenomicCodingEnd 
                     should be the last coding nucleotide's genomic coordinate 
                     on the exon. Otherwise, it should be the genomic 
                     coordinate of the first coding nucleotide.
StartPhase           Start phase of the exon. The position of an exon/intron 
                     boundary within a codon. A phase of zero means the 
                     boundary falls between codons, one means between the first 
                     and second base and two means between the second and third 
                     base. Exons have a start and end phase, whereas introns 
                     have just one phase. A boundary in a non-coding region has 
                     a phase of -1.
EndPhase             End phase of the exon.
NucleotideSequence   Nucleotide sequence for the exons. They can contain 
                     non-coding regions, e.g. UTRs. If the gene is in the 
                     reverse strand, the exon sequence should be the reverse 
                     complement of the genomic sequence.
==================== ===========================================================

You can find `user_transcript.csv`, an example of the required table in the 
`test/data` folder of the `thoraxe` repository at `GitHub`_. There is also 
this `Google spreadsheet`_ that you can modify and download later as a CSV file.


.. _GitHub: https://github.com/PhyloSofS-Team/thoraxe
.. _Google spreadsheet: https://docs.google.com/spreadsheets/d/1EEz1rsDCJdJeCl8jPoikTtskAKpEQiE3V8FGcCG55mk/edit?usp=sharing