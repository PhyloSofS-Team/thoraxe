Installation
============


Dependencies
------------


Clustal Omega
~~~~~~~~~~~~~

To run *ThorAxe* you need `Clustal Omega`_. Follow the instructions on that
link to install it.

If you are using **Ubuntu**, the easiest way to install *Clustal Omega* is:

::

   sudo apt-get install clustalo

If you are using **Windows 10**, you can use `clustalo` through *Ubuntu* with
the `Windows Subsystem for Linux`_.


NumPy
~~~~~

You can use the `pip`_ module of `Python 3`_ to install `NumPy`_ before
*ThorAxe*. The order is important because `numpy` should be installed before
`scikit-bio`_, a *ThorAxe* dependency:

::

   python -m pip install numpy

In *Ubuntu* you may need to use `python3` instead of `python`.

ThorAxe
-------

Finally, the easiest way to install *ThorAxe* is from PyPI_ using `pip`_ with
`Python 3`_:

::

   python -m pip install thoraxe


Installation from Source
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
the changes reflected into the installed version.


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
.. _Numpy: https://numpy.org/
.. _scikit-bio: http://scikit-bio.org/
.. _Clustal Omega: http://www.clustal.org/omega/
.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install-win10
.. _Docker image: https://hub.docker.com/r/diegozea/thoraxe
.. _Docker website: https://www.docker.com
.. _ipython3: https://ipython.readthedocs.io/en/stable/
.. _PyPI: https://pypi.org/project/thoraxe/
