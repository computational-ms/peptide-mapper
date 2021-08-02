#!/usr/bin/env python3
from setuptools import setup
import os

version_path = os.path.join(os.path.dirname(__file__), "peptide_mapper", "version.txt")
with open(version_path, "r") as version_file:
    mapper_version = version_file.read().strip()

req_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
with open(req_path, "r") as req_file:
    requirements = req_file.readlines()

setup(
    name="peptide_mapper",
    version=mapper_version,
    packages=["peptide_mapper"],
    package_dir={"peptide_mapper": "peptide_mapper"},
    package_data={
        "peptide_mapper": [
            "version.txt",
        ]
    },
    python_requires=">=3.7.0",
    install_requires=requirements,
    description="Peptide mapper",
    long_description="Mapper used to map peptides to a fasta file",
    author="",
    author_email="christian@fufezan.net",
    url="http://pymzml.github.com",
    license="The MIT license",
    platforms="any that supports python 3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: SunOS/Solaris",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
)
