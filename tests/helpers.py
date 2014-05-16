
import unittest
from mutgen import core
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def seqrecord(seqstring):
    "Silly little function for generating a simple seqrecord"
    sr = SeqRecord('whocares')
    sr.seq = Seq(seqstring)
    return sr


class TestGenerator(object):
    """Class for constructing TestCase subclasses with all the logic needed to test impirical frequencies
    based on a given mutation model."""

    def __init__(self, sample_size, almost_equal_places):
        """Construct the generator object; sample_size and almost_equal_places can be specified here and used
        as defaults for all TestCase subclasses generated"""

        self.sample_size = sample_size
        self.almost_equal_places = almost_equal_places

    def make_test(self, name, testseq, patterns, result_probs, sample_size=None, almost_equal_places=None):
        """Construct the new class, where `name` is the name of the new class, `testseq` is a sequence which
        will be tested, `patterns` specifies the mutation patterns (see `core.Mutator.__init__`), and
        `almost_equal_places` is a dictionary of kmers (or (kmer, new-base) pairs) to expected mutation
        frequencies."""

        def setUp(test_case):
            test_case.seqrecords = [seqrecord(testseq) for i in xrange(self.sample_size)]
            test_case.mutator = core.Mutator(patterns)
            test_case.sample_size = sample_size or self.sample_size
            test_case.almost_equal_places = almost_equal_places or self.almost_equal_places

        def test_frequencies(test_case):
            seqscan = test_case.mutator.seqscan(test_case.seqrecords)
            # Here key could either be a kmer or a (kmer, nt) pair
            results = dict((key, 0) for key in result_probs)
            for _, matching_kmers in seqscan:
                for result in matching_kmers:
                    kmer, mutated, mutated_to = map(lambda k: result[k], ['kmer', 'mutated', 'mutated_to'])
                    for key in [kmer, (kmer, mutated_to)]:
                        if mutated:
                            try:
                                results[key] += 1
                            except KeyError:
                                pass

            for kmer, count in results.iteritems():
                freq = float(count) / test_case.sample_size
                failmsg="Frequency of %s mutation was expected to be %s but was %s (diff: %s)." % (
                        kmer, result_probs[kmer], freq, abs(freq - result_probs[kmer]))
                test_case.assertAlmostEqual(freq, result_probs[kmer], places=test_case.almost_equal_places,
                        msg=failmsg)

        return type(name, (unittest.TestCase,), dict(setUp=setUp, test_frequencies=test_frequencies))


