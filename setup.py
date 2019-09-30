#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='thoraxe',
      description="Disentangle homology relationships between exons.",
      keywords=['homologous exons', 'protein isoforms', 'gene transcripts'],
      author="Diego Javier Zea, Hugues Richard & Elodie Laine",
      author_email='diegozea@gmail.com',
      url='https://github.com/PhyloSofS-Team/thoraxe',
      version='0.1.0',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English', 'Programming Language :: Python :: 3'
      ],
      entry_points={
          'console_scripts': [
              'transcript_query=thoraxe.transcript_query:main',
              'thoraxe=thoraxe:main'
          ],
      },
      packages=find_packages(),
      package_data={'thoraxe': ['subexons/assets/cluster_plots.html']},
      setup_requires=['pytest-runner'],
      install_requires=[
          'pandas', 'biopython>1.44', 'requests', 'numpy', 'recordclass',
          'scikit-bio', 'urllib3'
      ],
      test_suite='tests',
      tests_require=[
          'pytest', 'pytest-cov', 'coveralls', 'codecov', 'pytest-pylint',
          'pylint'
      ],
      include_package_data=True,
      zip_safe=False,
      license="MIT license")
