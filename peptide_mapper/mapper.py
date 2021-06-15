#!/usr/bin/env python
"""
usage:
    ./upeptide_mapper_1_0_0.py <input_file> <output_file>


"""

import sys
import os
import csv
import re
from collections import Counter, defaultdict
import itertools
from copy import deepcopy

try:
    import regex as regex

    finditer_kwargs = {"overlapped": True}
except:
    print("[ WARNING  ] Standard re module cannot find overlapping pattern")
    print("[   INFO   ] Consider installing the regex module")
    print("[   INFO   ] pip install -r requirements.txt")
    import re as regex

    finditer_kwargs = {}
import bisect

# increase the field size limit to avoid crash if protein merge tags
# become too long does not work under windows
if sys.platform != "win32":
    csv.field_size_limit(sys.maxsize)


def parse_fasta(io):
    """
    Small function to efficiently parse a file in fasta format.

    Keyword Arguments:
        io (obj): openend file obj (fasta file)

    Yields:
        tuple: fasta_id and sequence
    """
    id = None
    sequence = ""
    for line in io:
        line = line.strip()
        if line == "":
            continue
        if line[0] == ">":
            if id:
                yield (id, sequence)
            id = line[1:].strip()
            sequence = ""
        else:
            sequence += line
    if id:
        yield (id, sequence)


class UPeptideMapper:
    """
    UPeptideMapper V4

    Improved version of class version 3 (changes proposed by Christian)

    Note:
        Uses the implementation of Aho-Corasick algorithm pyahocorasick.
        Please refer to https://pypi.python.org/pypi/pyahocorasick/ for more
        information.



    """

    def __init__(self, fasta_database):

        # self.protein_indices            = {}
        self.protein_sequences = {}

        self.total_sequence_list = []

        self.protein_list = []
        self.protein_start_indices = []

        self.fasta_counter = 0
        self.len_total_sequence_string = 0

        self.peptide_2_protein_mappings = {}
        self.total_sequence_string = {}
        self.cache_database(fasta_database)

    def cache_database(self, fasta_database):
        """
        Function to cache the given fasta database.

        Args:
            fasta_database (str): path to the fasta database

        Note:

            If the same fasta_name is buffered again all info is purged from the
            class.
        """

        for protein_pos, (protein_id, seq) in enumerate(
            parse_fasta(open(fasta_database, "r").readlines())
        ):
            if protein_pos % 5000 == 0:
                print(
                    "[ upapa v4 ] Buffering protein #{0} of database {1}".format(
                        self.fasta_counter, os.path.basename(fasta_database)
                    ),
                    end="\r",
                )
            len_seq = len(seq)

            self.total_sequence_list.append(seq)

            self.protein_list.append(protein_id)
            self.protein_start_indices.append(
                self.len_total_sequence_string + protein_pos
            )  # protein start

            self.protein_sequences[protein_id] = seq
            self.len_total_sequence_string += len_seq  # offset for delimiter |
            self.fasta_counter += 1
        print()
        self.total_sequence_string = "|".join(self.total_sequence_list)
        return

    def map_peptides(self, peptide_list):
        """
        Function to map a given peptide list in one batch.

        Args:
            peptide_list (list): list with peptides to be mapped


        Returns:
            peptide_2_protein_mappings (dict): Dictionary containing
                peptides as keys and lists of protein mappings as values of the
                given fasta_name

        Note:
            Based on the number of peptides the returned mapping dictionary
            can become very large.

        Warning:
            The peptide to protein mapping is resetted if a new list o peptides
            is mapped to the same database (fasta_name).

        Examples::

            peptide_2_protein_mappings['PEPTIDE']  = [
                {
                    'start' : 1,
                    'end'   : 10,
                    'pre'   : 'K',
                    'post'  : 'D',
                    'id'    : 'BSA'
                }
            ]
        """
        import ahocorasick

        self.peptide_2_protein_mappings = defaultdict(list)

        # self.peptide_2_protein_mappings = defaultdict(list)
        self.automaton = ahocorasick.Automaton()
        for idx, peptide in enumerate(peptide_list):
            self.automaton.add_word(peptide, (idx, peptide))
        self.automaton.make_automaton()

        for match in self.automaton.iter(self.total_sequence_string):
            idx, (p_idx, m_peptide) = match
            len_m_peptide = len(m_peptide)

            protein_index = bisect.bisect(self.protein_start_indices, idx) - 1
            protein_name = self.protein_list[protein_index]
            protein_seq = self.protein_sequences[protein_name]
            start_in_total_seq_string = self.protein_start_indices[protein_index]
            stop_in_protein = idx - start_in_total_seq_string + 1
            start_in_protein = stop_in_protein - len(m_peptide)

            pre_pos = start_in_protein - 1

            if pre_pos < 0:
                pre = "-"
            else:
                pre = protein_seq[pre_pos]
            try:
                post = protein_seq[stop_in_protein]
            except:
                post = "-"

            self.peptide_2_protein_mappings[m_peptide].append(
                {
                    "start": start_in_protein + 1,
                    "end": stop_in_protein,
                    "pre": pre,
                    "post": post,
                    "id": protein_name,
                }
            )

        return self.peptide_2_protein_mappings


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    params = {
        "translations": {
            "modifications": [
                "M,opt,any,Oxidation",  # Met oxidation
                "C,fix,any,Carbamidomethyl",  # Carbamidomethylation
                "*,opt,Prot-N-term,Acetyl",  # N-Acteylation[]
            ],
            "aa_exception_dict": {
                "J": {
                    "original_aa": ["I", "L"],
                },
                "O": {
                    "original_aa": ["K"],
                    "unimod_name": "Methylpyrroline",
                },
            },
            "protein_delimiter": "<|>",
            "decoy_tag": "decoy_",
            "word_len": 6,
        },
        "prefix": None,
    }
    params["translations"]["database"] = sys.argv[3]
    params["translations"]["peptide_mapper_class_version"] = sys.argv[4]

    main(
        input_file=sys.argv[1],
        output_file=sys.argv[2],
        params=params,
    )
