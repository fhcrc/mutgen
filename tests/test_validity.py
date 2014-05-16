import unittest
from mutgen import core

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import random
random.seed("I can get you a toe dude. With nail polish!")

assert_almost_equal_places = 1
sample_size = 1000

def mut_pattern(motif, prob):
    return [motif, 0, [('A', prob), ('C', 0), ('G', 0), ('T', 0)]]


class TestSimpleSinglePattern(unittest.TestCase):
    def setUp(self):
        def seqrecord(seqstring):
            sr = SeqRecord('whocares')
            sr.seq = Seq(seqstring)
            return sr
        self.seqrecords = [seqrecord('CCCGTCCCCGACCCCCGGCCCCC') for i in xrange(sample_size)]
        self.mutator = core.Mutator(
                [mut_pattern('GT', 0.2),
                    mut_pattern('GA', 0.3),
                    mut_pattern('GG', 0.7)])

    def test_frequencies(self):
        seqscan = self.mutator.seqscan(self.seqrecords)
        results = dict(GT=0, GA=0, GG=0)
        for _, matching_kmers in seqscan:
            for kmer, mutated in matching_kmers:
                if mutated:
                    results[kmer] += 1

        self.assertAlmostEqual(float(results['GT']) / sample_size, 0.2, places=assert_almost_equal_places)
        self.assertAlmostEqual(float(results['GA']) / sample_size, 0.3, places=assert_almost_equal_places)
        self.assertAlmostEqual(float(results['GG']) / sample_size, 0.7, places=assert_almost_equal_places)


