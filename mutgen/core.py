"""
Core classes and control flow
"""

import csv
import copy
import re
import random
from Bio import SeqIO
from Bio.Seq import Seq


_AMBIGUOUS_MAP = {
       'R': 'GAR',      # GA
       'Y': 'TCY',      # TY
       'K': 'GTK',      # GT
       'M': 'ACM',      # AC
       'S': 'GCS',      # GC
       'W': 'ATW',      # AT
       'B': 'GTCYKSB',  # GTC (^A)
       'D': 'GATRKWD',  # GAT (^C)
       'H': 'ACTYMWH',  # ACT (^G)
       'V': 'GCARMSV',  # GCA (^T)
       'N': 'AGCTRYKMSWBDHVN'}


def regify_ambiguous(motif):
    def unambiguous(b):
        b = b.capitalize()
        amb = _AMBIGUOUS_MAP.get(b)
        return '[' + amb + ']' if amb else b
    return ''.join(map(unambiguous, motif))


class Mutif(object):
    def __init__(self, mutable_base, spec):
        """Takes an list of (base, prob) pairs, corresponding to the probability to mutating to a given base.
        The prob corresponding to `mutable_base` will be incremented by the value needed for the sum of the
        probs to equal 1."""

        self.mutable_base = mutable_base
        self.prob_ranges = []
        start = 0
        end = None
        prob_sum = sum(prob for _, prob in spec)
        assert prob_sum < 1.00000001, "Sum of probabilities must not be greater than 1.0"
        for base, prob in spec:
            if base == mutable_base:
                prob += 1.0 - prob_sum
            end = start + prob
            self.prob_ranges.append([start, end, base])
            start = end

    def __repr__(self):
        return "Mutif(%s)" % ", ".join("%s: %s" % (base, end - start)
                for start, end, base in self.prob_ranges)

    def mutate(self):
        p = random.random()
        _, _, bp = filter(lambda x: x[0] <= p and p < x[1], self.prob_ranges)[0]
        return bp


class KmerTooShort(Exception):
    """For some cleaner control flow excpetion handling"""
    pass


class Mutator(object):
    """Wrapper around a mutation model specification that takes care of the work of scanning sequences for
    mutations."""
    def __init__(self, spec):
        """`spec` arg should be a list of lists: `[motif, mut_index, [A, C, G, T]]`. Note that the mut
        index should only point to one of ACGT, and not to any ambiguous characters."""
        # Attributes:
        self.mut_index = max(i for _, i, _ in spec)
        self.width = self.mut_index + max(len(motif) - i for motif, i, _ in spec)
        self.mutable_bases = set()
        self.model = []

        # Normalize the motifs so the mut_index is the same for all
        def normalize(motif, mut_index):
            return '.' * (self.mut_index - mut_index) + motif
        spec = ((normalize(motif, mut_index), params) for motif, mut_index, params in spec)

        # Go through each motif and add to dictionary, but expand the motifs into
        for motif, params in spec:
            mutable_base = motif[self.mut_index]
            assert mutable_base.capitalize() in 'AGCT', "Mutable motif position must be non-ambiguous"
            self.mutable_bases.add(mutable_base)

            # Add our models
            self.model.append((regify_ambiguous(motif), Mutif(mutable_base, params)))

    @classmethod
    def from_file(cls, handle):
        """Load up mutator object from file."""
        reader = csv.DictReader(handle)
        return cls([[r['motif'], int(r['mut_index']), [(b, float(r.get(b, 0.0))) for b in 'ACGT']]
                for r in reader])

    def __repr__(self):
        return "Mutator(width: %s, mut_index: %s, components: %s, mutable_bases: %s)" % (self.width,
                self.mut_index, len(self.model), self.mutable_bases)

    def mutate(self, query_kmer):
        """Returns (mutated, mutated_to), where mutated is
          * None if bp corresponding to mut_index of query_kmer is not in the models mutable_bases; and
          * True/False otherwise, depending on whether there there is a mutation.
        Regardless of the value of mutated, mutated_to will be the new value of corresponding bp."""
        try:
            bp = query_kmer[self.mut_index]
        except IndexError:
            raise
        # First check to see if the bp is even mutable
        if bp in self.mutable_bases:
            # Check to see if any of the motifs match the kmer
            for motif, mutif in self.model:
                if re.match(motif, query_kmer):
                    # If they do, compute the new prob based on mutif
                    new_bp = mutif.mutate()
                    return (bp == new_bp, new_bp)
            # None of the motifs match, but since the bp is in mutable bases, we return False and the bp,
            # giving the interpretation that the model is assumed to be null here
            return (False, bp)
        else:
            # bp is not mutable so return (None, bp) as described in docstring
            return (None, bp)



def seqscan_iter(args):
    """Generator function that - given args.model, args.sequences and optionally args.width, yields a stream
    of pairs. Each pair contains the mutated sequence and a list of tuples that looks like (kmer, mutated),
    but only returns kmers where the mut_index position in the kmmer is in the models mutable_bases. The
    args.width specification controls the length of kmers that come out of this matching_kmers list.
    """
    mutator = Mutator.from_file(args.model)
    width = args.width or mutator.width
    seqrecords = SeqIO.parse(args.sequences, 'fasta')
    for sr in seqrecords:
        # Start new sequence with part not mutatable under model; we'll grow it
        seq = str(sr.seq)
        new_seq = seq[0:mutator.mut_index]
        # Init matching_mker list and go to town
        matching_kmers = []
        for i in xrange(mutator.mut_index, len(seq) - mutator.width + mutator.mut_index):
            kmer = seq[i:i+width]
            mutated, mutated_to = mutator.mutate(kmer)
            # this means the kmer's mutable bases is within the scope of our model, so we care about it for
            # the kmer out (matching_kmers)
            if mutated is not None:
                # Only add to list if kmer is actually as long as it's supposed to be
                if len(kmer) == width:
                    # XXX - need to hook in here to decide what to do with mutable base in output
                    #matching_kmers.append(dict(kmer=kmer, mutated=mutated))
                    matching_kmers.append((kmer, mutated))
            # Add the mutated_to value (which may just actually be the old bp) to the growing seq
            new_seq += mutated_to
        # add the last bit of the sequence on, the part on the other side of mut_index that can't be mutated
        new_seq += seq[len(seq) - mutator.width + mutator.mut_index:]
        # Create a new seqrecord copy and give it the new_seq
        new_sr = copy.copy(sr)
        new_sr.seq = Seq(new_seq)
        # Yield the mutated seq record and the matching kmers
        yield new_sr, matching_kmers

 
def seqgen_handler(args):
    raise ValueError("The seqgen functionality hasn't actually been built out yet. Sorry for the tease.")



