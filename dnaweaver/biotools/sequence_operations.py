from io import StringIO
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.
    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"CGGCAT"``.
    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def random_dna_sequence(length, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility.
    """
    if seed is not None:
        np.random.seed(seed)
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def load_record(filename, name="unnamed"):
    """Load a sequence file (genbank or fasta) as a Biopython record."""
    if filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError("Unknown format for file: %s" % filename)
    record.id = name
    record.name = name.replace(" ", "_")[:20]
    return record


def string_to_sequence(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    if (matches is not None) and (matches.groups()[0] == string):
        return (string, "ATGC")

    for fmt in ("fasta", "genbank"):
        try:
            stringio = StringIO(string)
            return (str(SeqIO.read(stringio, fmt).seq), fmt)
        except Exception:
            pass
    raise ValueError("Invalid sequence format")


def file_to_sequence(filename):
    """Import a file in fasta, genbank... as a simple ATGC string."""
    with open(filename, "r") as f:
        return string_to_sequence(f.read())


def sequence_to_atgc(seq):
    if isinstance(seq, str):
        return seq
    if hasattr(seq, "seq"):
        return str(seq.seq)
    return str(seq)
