Quick Introduction
==================

*ThorAxe*  finds *s-exons*, i.e. potential orthologous exonic regions, and 
allows the creation of *Evolutionary Splicing Graphs (ESG)*. There are more 
details in our `publication`_. It offers three programs towards this goal: 
`transcript_query`, `add_transcripts` and `thoraxe`. The three programs are 
meant to be run in that order, being `add_transcripts` optional:

.. tip::
    You can check our `YouTube tutorial`_ to see how you can use *ThorAxe* and 
    visualise the results.

I. Download transcript information from Ensembl
-----------------------------------------------

`transcript_query` downloads the information that *ThorAxe* needs from Ensembl; 
it only requires a gene name. You can find more information on our 
documentation for `transcript_query`_ or by running `transcript_query -h`. 
For example, let's say you want to analyse *MAPK8* between human and mouse:

::

    transcript_query MAPK8 --specieslist homo_sapiens,mus_musculus

[OPTIONAL] Add transcripts to the previously download data
-----------------------------------------------------------

If you have transcripts that are not available at Ensembl, e.g. data coming 
from your experiments, you can manually add them using `add_transcripts`_.

II. Run ThorAxe
----------------

Finally, you need to run `thoraxe` on that data. There are multiple 
configuration options that you can check in the `thoraxe`_ documentation or by 
running `thoraxe -h`.

::

    thoraxe -i MAPK8


.. _publication: https://doi.org/10.1101/2020.11.14.382820
.. _transcript_query: https://phylosofs-team.github.io/thoraxe/programs/transcript_query.html
.. _add_transcripts: https://phylosofs-team.github.io/thoraxe/programs/add_transcripts.html
.. _thoraxe: https://phylosofs-team.github.io/thoraxe/programs/thoraxe.html
.. _YouTube tutorial: https://www.youtube.com/watch?v=Z96985kX-uY