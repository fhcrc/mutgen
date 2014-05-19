
# MutGen - motif based mutation simulation library

MutGen allows you to mutate sequences according to a specified mutation model.
You can currently scan a set of sequences, generating mutations as you go according to the model (see the `seqscan` cli subcommand).
We hope to soon also implement generation of random kmers/sequences, which get similarly mutated according to the model.

## Model files

To use MutGen you must specify a motif-based mutation model.
This can be done via a CSV file that contains the columns `motif`, `mut_index` and one or more of `"ACGT"`.

The given `motif` can contain ambiguous characters in any position not corresponding to the `mut_index` (the _mutable base_; that is, the position in the motif that is mutatable).
The `ACGT` columns represent the probability for a given motif that the mutable base will mutate to the corresponding base.

If, for example, the mutable base for some motif is `G`, the probability specified in a `G` column will be ignored - the probability of leaving not mutating the base is taken as `1 - sum(other_probs)`.
If there is not a column for a given base, the probability of mutating to that base is assumed to be 0.
Motifs not represented in the specification are assumed to have 0 probability of mutation.

For an example model file, please see `examples/model.csv`

## CLI

For a complete set of features type `mutgen -h`, or `mutgen <subcommand> -h` for the `seqscan` or `seqgen` subcommands.

### OUTPUTS

The seqscan function has two optional output files.
The `-s|--out-seqs` file will contain the mutated sequences (with mutations "in place").
The `-k|--out-kmers` file will contain a list of the kmers with a mutable base position matching that of the model specification, and a boolean column indicating whether the kmer was found to be mutated under the model.

## Deep Thoughts...

### Future:

* add cols for mutable_base, new_base, and new_kmer; swap out mutated_to for new_base
* add option for fully expanding specified patterns into dicts or not
* option for specifying motifs directly using [] regexps?

### Valiate?:

* no overlapping motifs - currently, the software just matches on the first pattern and mutates that.
  Should either let it match multiple and add the probs or maybe just validate that none of the motifs overlap (though this might be rather tricky)

