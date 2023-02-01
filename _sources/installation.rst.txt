Installation
============


Dependencies
------------

*ThorAxe* needs **Python version 3.6.1 or higher, but less than 3.9** and:

ProGraphMSA
~~~~~~~~~~~

*ThorAxe* needs `ProGraphMSA`_ to perform the multiple sequence alignments. 
Follow the instructions on that link to install it or download one of the 
following compiled binaries:

- `ProGraphMSA Linux 32 bits`_
- `ProGraphMSA Linux 64 bits`_
- `ProGraphMSA OSX`_

It is important to rename the program to `ProGraphMSA` and to put it in the 
path. If you are using **Windows 10**, you can use `ProGraphMSA` through 
*Ubuntu* with the `Windows Subsystem for Linux`_. For example, to install the 
64 bits version in *Ubuntu* you can do:

::

   wget  https://github.com/PhyloSofS-Team/ProGraphMSA/raw/master/bin/ProGraphMSA_64
   chmod a+x ProGraphMSA_64
   sudo mv ProGraphMSA_64 /bin/ProGraphMSA


ThorAxe
-------

Finally, the easiest way to install *ThorAxe* is from PyPI_ using `pip`_ with
`Python 3`_:

::

   python -m pip install thoraxe

In *Ubuntu*, you may need to use `python3` instead of `python`.

::

   python3 -m pip install thoraxe

Installation from source
------------------------

If you want to install the latest and in-development version of *ThorAxe*, clone
its source code from the `thoraxe` `GitHub repo`_ using `git`_:

::

   git clone https://github.com/PhyloSofS-Team/thoraxe.git

Then, you can use the `pip`_ module of `Python 3`_ to install it:

::

   cd thoraxe
   python -m pip install -e .

The `-e` flag allows you to update the `thoraxe` repo using `git pull` and get
the changes reflected in the installed version.


Docker image
------------

We offer a `Docker image`_ with *ThorAxe* and all its dependencies installed.
To run this image, you need to install *Docker* from the `Docker website`_.

The following example is going to run the *ThorAxe Docker* image using
*Windows PowerShell*. The local directory is going to be bind mounted in
`/project`. The actual folder is `${PWD}` in *Windows PowerShell* and `$(pwd)`
in *Unix*. Please note that in *Ubuntu*, you may need `sudo` to run `docker`.

::

   docker run -ti --rm --mount type=bind,source=${PWD},target=/project diegozea/thoraxe


After this, we have access to the `bash` terminal of an *Ubuntu 18.04* image
with *ThorAxe*, `tree` to help explore the results and `ipython3`_ to make
easier the use of the `thoraxe` *Python* module.

.. figure :: _static/docker.gif

   ThorAxe example using Docker in Windows 10.


.. _git: https://git-scm.com/
.. _GitHub repo: https://github.com/PhyloSofS-Team/thoraxe
.. _pip: https://pip.pypa.io/en/stable/installing/
.. _Python 3: https://www.python.org/
.. _ProGraphMSA: https://github.com/acg-team/ProGraphMSA
.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install-win10
.. _Docker image: https://hub.docker.com/r/diegozea/thoraxe
.. _Docker website: https://www.docker.com
.. _ipython3: https://ipython.readthedocs.io/en/stable/
.. _PyPI: https://pypi.org/project/thoraxe/
.. _ProGraphMSA Linux 32 bits: https://github.com/PhyloSofS-Team/ProGraphMSA/raw/master/bin/ProGraphMSA_32
.. _ProGraphMSA Linux 64 bits: https://github.com/PhyloSofS-Team/ProGraphMSA/raw/master/bin/ProGraphMSA_64
.. _ProGraphMSA OSX: https://github.com/PhyloSofS-Team/ProGraphMSA/raw/master/bin/ProGraphMSA_osx