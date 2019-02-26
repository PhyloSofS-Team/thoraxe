#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='exonhomology',
    description="Disentangle homology relationships between exons.",
    keywords=['homologous exons', 'protein isoforms', 'gene transcripts'],
    author="Diego Javier Zea, Hugues Richard & Elodie Laine",
    author_email='diegozea@gmail.com',
    url='https://github.com/PhyloSofS-Team/exonhomology',
    version='0.1.0',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points={
        'console_scripts': [
            'transcript_query=exonhomology.transcript_query:main',
            'exonhomology=exonhomology:main'
        ],
    },
    packages=find_packages(include=['exonhomology']),
    package_data={'exonhomology': ['subexons/assets/cluster_plots.html']},
    setup_requires=['pytest-runner'],
    install_requires=[
        'pandas', 'biopython', 'requests', 'numpy', 'recordclass', 'scikit-bio'
    ],
    test_suite='tests',
    tests_require=[
        'pytest', 'pytest-cov', 'coveralls', 'codecov', 'pytest-pylint'
    ],
    include_package_data=True,
    zip_safe=False,
    license="MIT license")
