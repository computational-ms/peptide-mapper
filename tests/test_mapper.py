#!/usr/bin/env python3
# encoding: utf-8

import unittest
import pprint
import sys
from pathlib import Path

# from ursgal import umapmaster as umama
import pprint
import os
from peptide_mapper.mapper import UPeptideMapper

try:
    import regex

    regex_module_installed = True
except:
    regex_module_installed = False


TEST_FASTA = [
    ">Protein1\n",
    "ELVISLIVES\n",
    "\n",
    ">Protein2\n",
    "KLEINERPENNER\n",
    ">Protein3\n",
    "WHYISELVIS\n",
    "HELEAVING\n",
    ">T1",
    "AAAAAAA___AAAAAAAAAA_____AAAAA__________AAAAAAAAAA",
    # ----+----|----+----|----+----|----+----|----+----|
    ">Overlapping",
    "GGGGGGGGGG",
    ">GnomeChompsky",
    "MUSTACHIO",
    ">Overlap_protein",
    "OVERLAP",
    ">Overlap_protein2",
    "OVERLAP",
]

TEST_FASTA_TWO = [
    ">Ze_other_protein2\n",
    "KLEINERRETTER\n",
]


class UMapMaster(unittest.TestCase):
    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def setUp(self):
        self.database_path = Path(__file__).parent / "data" / "Test.fasta"
        with open(self.database_path, "w") as io:
            for line in TEST_FASTA:
                print(line.strip(), file=io)
        self.upapa_class = UPeptideMapper(fasta_database=self.database_path)
        return

    # @unittest.skipIf(
    #     sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    # )
    # def test_purge_fasta(self):
    #     self.upapa_class.purge_fasta_info("Test.fasta")
    #     self.assertFalse("Test.fasta" in self.upapa_class.protein_sequences.keys())
    #     self.assertFalse("Test.fasta" in self.upapa_class.protein_indices.keys())
    #     self.assertFalse(
    #         "Test.fasta" in self.upapa_class.peptide_2_protein_mappings.keys()
    #     )
    #     return

    # @unittest.skipIf(
    #     sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    # )
    # def test_fasta_id_parsed_and_available(self):
    #     input_fastas = []
    #     for line in TEST_FASTA:
    #         if line.startswith(">"):
    #             input_fastas.append(line[1:].strip())

    #     self.assertEqual(
    #         sorted(self.upapa_class.protein_sequences["Test.fasta"].keys()),
    #         sorted(input_fastas),
    #     )
    #     return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_peptide_mapping_1(self):
        self.upapa_class.map_peptides(["ELVIS"])
        self.assertEqual(
            [
                (d["id"], d["start"])
                for d in self.upapa_class.peptide_2_protein_mappings["Test.fasta"][
                    "ELVIS"
                ]
            ],
            [
                ("Protein1", 1),
                ("Protein3", 6),
            ],
        )

        self.upapa_class.map_peptides(["MUSTACHIO"])
        self.assertEqual(
            [
                (d["id"], d["start"])
                for d in self.upapa_class.peptide_2_protein_mappings["Test.fasta"][
                    "MUSTACHIO"
                ]
            ],
            [
                ("GnomeChompsky", 1),
            ],
        )
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_peptide_mapping_2(self):
        maps = self.upapa_class.map_peptides(["VISHE"])["VISHE"]
        self.assertEqual(len(maps), 1)
        for mapping in maps:
            self.assertEqual("Protein3", mapping["id"])
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_incremental_cache_buildups(self):
        """"""
        self.upapa_class.cache_database(self.database_path)
        # map with one parsed fasta
        maps = self.upapa_class.map_peptides(["KLEINER"])["KLEINER"]
        breakpoint()
        self.assertEqual(len(maps), 1)
        tmp_database_path = os.path.join("tests", "data", "tmp_Test.fasta")
        with open(tmp_database_path, "w") as io:
            for line in TEST_FASTA_TWO:
                print(line.strip(), file=io)

        self.upapa_class.cache_database(tmp_database_path, "tmp_Test.fasta")
        maps = self.upapa_class.map_peptides(["KLEINER"])["KLEINER"]
        breakpoint()
        self.assertEqual(len(maps), 1)
        self.assertEqual(
            len(
                self.upapa_class.map_peptides(["KLEINER"], "tmp_Test.fasta")["KLEINER"]
            ),
            1,
        )

        os.remove(tmp_database_path)
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_peptide_gt_word_but_not_continous(self):
        self.assertEqual(self.upapa_class.map_peptides(["WHYELVIS"])["WHYELVIS"], [])
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_multiple_occurrence_in_one_seq(self):
        maps = self.upapa_class.map_peptides(["AAAAAAAAAA"])["AAAAAAAAAA"]
        self.assertEqual(len(maps), 2)
        sorted_maps = sorted(maps, key=lambda x: x["start"])
        self.assertEqual(
            sorted(maps, key=lambda x: x["start"]),
            sorted(
                [
                    {"pre": "_", "post": "_", "start": 11, "end": 20, "id": "T1"},
                    {"pre": "_", "post": "-", "start": 41, "end": 50, "id": "T1"},
                ],
                key=lambda x: x["start"],
            ),
        )
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_multiple_occurrence_with_opverlap_in_one_seq(self):
        maps = self.upapa_class.map_peptides(["GGGGGGG"])["GGGGGGG"]
        # print(maps)
        self.assertEqual(len(maps), 4)
        self.assertEqual(sorted([m["start"] for m in maps]), [1, 2, 3, 4])
        maps = self.upapa_class.map_peptides(["GGGG"])["GGGG"]
        # print(maps)
        self.assertEqual(len(maps), 7)
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_overlap_produces_no_match(self):
        maps = self.upapa_class.map_peptides(["LAPOVER"])["LAPOVER"]
        self.assertEqual(len(maps), 0)
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_start_and_end_of_sequence(self):
        maps = self.upapa_class.map_peptides(["OVERLAP"])["OVERLAP"]
        # print(maps)
        for map_dict in maps:
            self.assertEqual(map_dict["post"], "-")
            self.assertEqual(map_dict["pre"], "-")
        return

    @unittest.skipIf(
        sys.platform == "win32", "pyahocorasick not installed via pip in Windows"
    )
    def test_end_of_sequence_and_fasta(self):
        maps = self.upapa_class.map_peptides(["LAP"])["LAP"]
        # print(maps)
        for map_dict in maps:
            self.assertEqual(map_dict["post"], "-")
        return


if __name__ == "__main__":
    unittest.main()
