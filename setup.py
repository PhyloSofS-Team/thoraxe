#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    author="Diego Javier Zea, Hugues Richard & Elodie Laine",
    author_email='diegozea@gmail.com',
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
    description="Disentangle homology relationships between exons.",
    entry_points={
        'console_scripts': [
            'transcript_query=exonhomology.transcript_query:main',
            'exonhomology=exonhomology:main'
        ],
    },
    install_requires=[
        'numpy', 'pandas', 'matplotlib', 'seaborn', 'biopython', 'requests',
        'scikit-bio', 'distro', 'recordclass'
    ],
    license="MIT license",
    include_package_data=True,
    keywords=['homologous exons', 'protein isoforms', 'gene transcripts'],
    name='exonhomology',
    packages=find_packages(include=['exonhomology']),
    setup_requires=['pytest-runner'],
    test_suite='tests',
    tests_require=[
        'pytest', 'pytest-cov', 'coveralls', 'codecov', 'pytest-pylint'
    ],
    url='https://github.com/PhyloSofS-Team/exonhomology',
    version='0.1.0',
    zip_safe=False,
)
