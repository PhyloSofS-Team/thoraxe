transcript_query
================

.. argparse::
   :ref: thoraxe.transcript_query.transcript_query.parse_command_line
   :prog: transcript_query


.. caution::
   It is possible to get more transcripts and genes by allowing `1:n` or `n:m` 
   orthology using the `--orthology` argument. However, `thoraxe` can not 
   ensure that the *s-exons* contain orthologous exonic regions under 
   that setup.