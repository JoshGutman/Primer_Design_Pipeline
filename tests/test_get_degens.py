import unittest
import sys
import os

from get_degens import *
from get_primers import Primer

class TestReverseComplement(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("T"), "A")
        self.assertEqual(reverse_complement("C"), "G")
        self.assertEqual(reverse_complement("G"), "C")

    def test_sequence(self):
        self.assertEqual(reverse_complement("ATAGGTACC"), "GGTACCTAT")

    def test_lowercase_sequence(self):
        self.assertEqual(reverse_complement("acgt"), "ACGT")

    def test_degens(self):
        self.assertEqual(reverse_complement("RYSWKMBVDHN"), "NDHBVKMWSRY")

    def test_invalid_char(self):
        with self.assertRaises(ValueError):
            reverse_complement("ACG?")

    def test_empty_seq(self):
        self.assertEqual(reverse_complement(""), "")


class TestNumCombos(unittest.TestCase):

    def test_num_combos(self):
        self.assertEqual(num_combos("ACGT"), 1)
        self.assertEqual(num_combos("ACGM"), 2)
        self.assertEqual(num_combos("ACGV"), 3)
        self.assertEqual(num_combos("ACGN"), 4)
        self.assertEqual(num_combos("MRWS"), 16)
        self.assertEqual(num_combos("MRWV"), 24)
        self.assertEqual(num_combos("MRWN"), 32)

    def test_lowercase_seq(self):
        self.assertEqual(num_combos("acgt"), 1)

    def test_invalid_char(self):
        with self.assertRaises(ValueError):
            num_combos("ACG?")

    def test_empty_seq(self):
        self.assertEqual(num_combos(""), 1)


class TestAmbiguityCode(unittest.TestCase):

    def test_ambiguity_code(self):
        self.assertEqual(ambiguity_code(["A"]), "A")
        self.assertEqual(ambiguity_code(["A", "G"]), "R")
        self.assertEqual(ambiguity_code(["C", "T"]), "Y")
        self.assertEqual(ambiguity_code(["C", "G"]), "S")
        self.assertEqual(ambiguity_code(["A", "T"]), "W")
        self.assertEqual(ambiguity_code(["G", "T"]), "K")
        self.assertEqual(ambiguity_code(["A", "C"]), "M")
        self.assertEqual(ambiguity_code(["C", "G", "T"]), "B")
        self.assertEqual(ambiguity_code(["A", "G", "T"]), "D")
        self.assertEqual(ambiguity_code(["A", "C", "T"]), "H")
        self.assertEqual(ambiguity_code(["A", "C", "G"]), "V")
        self.assertEqual(ambiguity_code(["A", "C", "G", "T"]), "N")

    def test_unsorted_seq(self):
        self.assertEqual(ambiguity_code(["T", "G", "C", "A"]), "N")

    def test_lowercase_seq(self):
        self.assertEqual(ambiguity_code(["a", "c", "g", "t"]), "N")

    def test_invalid_char(self):
        with self.assertRaises(ValueError):
            ambiguity_code(["A", "C", "G", "?"])

    def test_empty_seq(self):
        self.assertEqual(ambiguity_code([]), "N")


class TestGetDegens(unittest.TestCase):

    def test_get_degens(self):
        genomes = [
            "ACGTACGTAC",
            "CCGTACGTAC",
            "GCGTACGTAC",
            "TCGTACGTAC",
        ]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, 0)

        self.assertEqual(forward.sequence, "NCGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGN")

    def test_ignore_percent(self):
        genomes = [
            "ACGTACGTAC",
            "ACGTACGTAC",
            "CCGTACGTAC",
            "GCGTACGTAC",
            "TCGTACGTAC",
        ]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, .4)
        
        self.assertEqual(forward.sequence, "ACGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGT")

    def test_ignore_percent_greater_than_1(self):
        genomes = [
            "ACGTACGTAC",
            "ACGTACGTAC",
            "CCGTACGTAC",
            "GCGTACGTAC",
            "TCGTACGTAC",
        ]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, 40)
        
        self.assertEqual(forward.sequence, "ACGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGT")

    def test_primer_sequence_not_found(self):
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        
        genomes = ["----------"]
        primers = [Primer("LEFT", "ACGTACGTAC", 0), Primer("RIGHT", "GTACGTACGT", 0)]
        get_degens(primers, genomes, 0)
        
        sys.stdout.close()
        sys.stdout = old_stdout
        
        self.assertEqual(len(primers), 0)

    def test_primer_sequence_too_degenerate(self):
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        
        genomes = [
            "ACGTACGTAC",
            "ACGTA-----",
        ]
        primers = [Primer("LEFT", "ACGTACGTAC", 0), Primer("RIGHT", "GTACGTACGT", 0)]
        get_degens(primers, genomes, .99)
        
        sys.stdout.close()
        sys.stdout = old_stdout
        
        self.assertEqual(len(primers), 0)

    def test_lowercase_seq_in_primer(self):
        genomes = ["ACGTACGTAC"]
        forward = Primer("LEFT", "acgtacgtac", 0)
        reverse = Primer("RIGHT", "gtacgtacgt", 0)
        get_degens([forward, reverse], genomes, 0)
        
        self.assertEqual(forward.sequence, "ACGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGT")

    def test_invalid_char_in_primer(self):
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        
        genomes = ["ACGTACGTAC"]
        forward = Primer("LEFT", "?CGTACGTAC", 0)
        reverse = Primer("RIGHT", "?TACGTACGT", 0)        
        
        with self.assertRaises(ValueError):
            get_degens([forward, reverse], genomes, 0)

        sys.stdout.close()
        sys.stdout = old_stdout

    def test_degen_in_genome(self):
        genomes = [
            "ACGTACGTAC",
            "NCGTACGTAC",
            "-CGTACGTAC",
        ]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, .5)
        
        self.assertEqual(forward.sequence, "NCGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGN")
        
    def test_invalid_char_in_genome(self):
        genomes = [
            "ACGTACGTAC",
            "?CGTACGTAC",
        ]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, 0)
        
        self.assertEqual(forward.sequence, "ACGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGT")

    def test_lowercase_seq_in_genome(self):
        genomes = ["acgtacgtac"]
        forward = Primer("LEFT", "ACGTACGTAC", 0)
        reverse = Primer("RIGHT", "GTACGTACGT", 0)
        get_degens([forward, reverse], genomes, 0)
        
        self.assertEqual(forward.sequence, "ACGTACGTAC")
        self.assertEqual(reverse.sequence, "GTACGTACGT")


if __name__ == "__main__":
    unittest.main()
