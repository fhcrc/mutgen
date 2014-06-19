#!/usr/bin/env python
"Command line interface to mutgen"

import argparse
import random
import csv              
from core import Mutator
from Bio import SeqIO


class PositivesReached(Exception):
    pass


def seqscan_handler(args):
    mutator = Mutator.from_file(args.model)
    if args.width:
        width = args.width
    elif args.extra:
        width = mutator.width + args.extra
    else:
        # That is, default to mutator.width
        width = None

    positives = 0

    kmer_fieldnames = ['kmer', 'mutated', 'mutated_to', 'sequence']
    if args.replicates > 1:
        kmer_fieldnames.append('replicate')
    if args.positives:
        kmer_fieldnames.append('positives')

    kmer_writer = csv.DictWriter(args.out_kmers, fieldnames=kmer_fieldnames, extrasaction='ignore') if args.out_kmers else None
    kmer_writer.writeheader()

    try:
        for replicate in xrange(args.replicates):
            args.sequences.seek(0)
            sequences = SeqIO.parse(args.sequences, 'fasta')

            for mutated_seq, matching_kmers in mutator.seqscan(sequences, width=width,
                    match_mutated=not args.relaxed_sampling):
                # First write out sequence if needed
                if args.out_seqs:
                    SeqIO.write([mutated_seq], args.out_seqs, 'fasta')
                # Next write out kmers if needed
                if args.out_kmers:
                    for kmer_data in matching_kmers:
                        mutated = kmer_data['mutated']
                        if mutated:
                            positives += 1
                        # Break out if we have a max_positives set
                        if args.positives and positives > args.positives:
                            raise PositivesReached
                        # If we haven't reached positives_max yet, write row to file
                        kmer_data.update(positives=args.positives, replicate=replicate, sequence=mutated_seq.id)
                        kmer_writer.writerow(kmer_data)

    # Break out if we are stopping at a given positives count
    except PositivesReached:
        pass

    # Close up shop
    finally:
        if args.out_seqs:
            args.out_seqs.close()
        if args.out_kmers:
            args.out_kmers.close()


def seqgen_handler(args):
    raise ValueError("The seqgen functionality hasn't actually been built out yet. Sorry for the tease.")


def setup_common_args(subparser):
    subparser.add_argument('model', type=argparse.FileType('r'),
        help="""Input CSV with columns motif,mut_index,A,C,G,T specifying the probability that the mut_index of
        motif will mutate to the given base pair. A,C,G,T columns left blank imply zero probability. The
        probability of mutating a base pair to itself is ignored, and set to 1 - sum(other_probs).
        The base corresponding to mut_index must be a non-ambiguous character for all motifs.""")
    subparser.add_argument('-w', '--width', type=int, help="""Specify width of generated k-mer in seqgen, or
        in seqscan, the length of the k-mers considered for mutation and written to --motifs output file,
        if specified.""")
    subparser.add_argument('-e', '--extra', type=int, help="""Specify width by specifying how much extra
            padding should be added to the natural width of the model specification.""")
    subparser.add_argument('-r', '--replicates', type=int, default=1,
        help="""Repeat entire process the specified number of times.""")
    subparser.add_argument('-S', '--seed', help="Specify random seed")


def setup_seqscan(subparsers):
    subparser = subparsers.add_parser('seqscan', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_common_args(subparser)
    subparser.add_argument('sequences', type=argparse.FileType('r'), help="Input Fasta with sequences to mutate")
    subparser.add_argument('-s', '--out-seqs', type=argparse.FileType('w'),
        help="Output mutated sequences")
    subparser.add_argument('-k', '--out-kmers', type=argparse.FileType('w'),
        help="Output CSV of mutated and non-mutated kmers")
    subparser.add_argument('-p', '--positives', type=int,
        help="""If set, scan through sequences till specified number of positives have been hit (does not
        currently 'rescan' sequences, but may in future)""")
    subparser.add_argument('-R', '--relaxed-sampling', action="store_true",
        help="""By default, this program 'scans' along the sequence, mutating things as it goes. This means
        that the motif surrounding a mutable base depends on the mutations that came before it. Specifying
        this flag relaxes this assumption, so that the pre-mutation context is all that's considered in
        computating mutation probabilities.""")
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


def main():
    args = get_args()
    if args.seed:
        random.seed(args.seed)
    args.func(args)


if __name__ == '__main__':
    main()


