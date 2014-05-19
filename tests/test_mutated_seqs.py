
import unittest
import random
import helpers

from mutgen import core

random.seed("Say what you want about the tenets of national socialism dude; at least it's an ethos!")


class TestZeroMutPosition(unittest.TestCase):
    def setUp(self):
        self.model = core.Mutator([('TGA', 0, {'A': 0.99999})])
        self.seqs = [helpers.seqrecord("TTTTGATTTT")]

    def test_mutated_seq(self):
        for mutated_seq, _ in self.model.seqscan(self.seqs):
            self.assertEqual(str(mutated_seq.seq), "TTTAGATTTT")


class TestNonzeroPosition(unittest.TestCase):
    def setUp(self):
        self.model = core.Mutator([('TGA', 1, {'A': 0.99999})])
        self.seqs = [helpers.seqrecord("TTTTGATTTT")]

    def test_mutated_seq(self):
        for mutated_seq, _ in self.model.seqscan(self.seqs):
            self.assertEqual(str(mutated_seq.seq), "TTTTAATTTT")

class TestCustomWidth(unittest.TestCase):
    def setUp(self):
        self.model = core.Mutator([('TGA', 1, {'A': 0.99999})])
        self.seqs = [helpers.seqrecord("TTTTGATTTT")]

    def test_mutated_seq(self):
        for mutated_seq, _ in self.model.seqscan(self.seqs, width=4):
            self.assertEqual(str(mutated_seq.seq), "TTTTAATTTT")

