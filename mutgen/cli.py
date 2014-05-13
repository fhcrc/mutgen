#!/usr/bin/env python
"Command line interface to mutgen"

import argparse
import random
import csv
import itertools as it
from core import Mutator, seqscan_iter, seqgen_handler
from utils import iunzip
from Bio import SeqIO


def seqscan_handler(args):
    mutated_seqs, matching_kmers = iunzip(seqscan_iter(args))
    # First write out sequences if needed
    if args.out_seqs:
        SeqIO.write(mutated_seqs, args.out_seqs, 'fasta')
    # Next write out kmers if needed
    if args.out_kmers:
        # XXX - we should throw in the sequence name as well
        #writer = csv.DictWriter(args.out_kmers)
        writer = csv.writer(args.out_kmers)
        writer.writerow(['kmer', 'mutated'])
        writer.writerows(it.chain(*matching_kmers))


def setup_common_args(subparser):
    subparser.add_argument('model', type=argparse.FileType('r'),
        help="""Input CSV with columns motif,mut_index,A,C,G,T specifying the probability that the mut_index of
        motif will mutate to the given base pair. A,C,G,T columns left blank imply zero probability. The
        probability of mutating a base pair to itself is ignored, and set to 1 - sum(other_probs).
        The base corresponding to mut_index must be a non-ambiguous character for all motifs.""")
    subparser.add_argument('-w', '--width', type=int, help="""Specify width of generated k-mer in seqgen, or
        in seqscan, the length of the k-mers considered for mutation and written to --motifs output file,
        if specified.""")
    subparser.add_argument('-r', '--replicates', type=int,
        help="""Repeat entire process the specified number of times.""")
    subparser.add_argument('-S', '--seed', help="Specify random seed")


def setup_seqscan(subparsers):
    subparser = subparsers.add_parser('seqscan', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_common_args(subparser)
    subparser.add_argument('sequences', type=argparse.FileType('r'), help="Input Fasta with sequences to mutate")
    subparser.add_argument('-s', '--out-seqs', type=argparse.FileType('w'),
        help="Output sequences with inline mutations")
    subparser.add_argument('-k', '--out-kmers', type=argparse.FileType('w'),
        help="Output CSV of mutated and non-mutated kmers")
    subparser.add_argument('-p', '--positives', type=int,
        help="""If set, scan through sequences till the specified number of positives have been hit""")
    subparser.set_defaults(func=seqscan_handler)


def setup_seqgen(subparsers):
    subparser = subparsers.add_parser('seqgen', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_common_args(subparser)
    subparser.add_argument('output', type=argparse.FileType('w'),
        help="""Output CSV of mutated and non-mutated motifs""")
    subparser.add_argument('-p', '--positives', type=int,
        help="""Generate motifs and mutate until the specified number of positions have been mutated""")
    subparser.add_argument('-b', '--background-model',
        help="""Pymix model file for simulating sequences for mutation consideration""")
    subparser.set_defaults(func=seqgen_handler)


def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title="subcommands",
        help="""Select operation mode (for help with specified mode, add -h flag)""")
    setup_seqscan(subparsers)
    setup_seqgen(subparsers)
    return parser.parse_args()


def main(args):
    if args.seed:
        random.seed(args.seed)
    args.func(args)


def test(args):
    m = Mutator.from_file(args.model)
    print m
    print m.model
    print m.mutate('CAAGCTK')


if __name__ == '__main__':
    main(get_args())


