"""Classes of constraints to model cloning capabilities in
CommercialDnaOffer, DnaAssemblyMethod, etc.
"""

from Bio import Restriction
import re
from ..biotools import gc_content, reverse_complement


class NoPatternConstraint:
    """Constraint class forbidding a given pattern in DNA sequences.

    Class of callables (sequence)-> True/False whether the sequence contains
    the pattern.

    Can be useful for defining constraints in DNA assembly methods or
    DNA providers.

    The interest of having this as a class is that a DnaSupplier using this
    constraint can be displayed as a string with the pattern appearing
    explicitly, which would not be the case for a function

    Parameters
    ----------

    pattern=None, enzyme=None, is_regex=False, with_revcomp=True
    """

    def __init__(self, pattern=None, enzyme=None, is_regex=False, with_revcomp=True):

        self.biopython_enzyme = None
        if enzyme is not None:
            if enzyme in Restriction.__dict__:
                biopython_enzyme = Restriction.__dict__[enzyme]
                if all([c in "ATGC" for c in biopython_enzyme.site]):
                    pattern = biopython_enzyme.site
                else:
                    self.biopython_enzyme = biopython_enzyme
            else:
                raise ValueError("Unknown enzyme: %s" % enzyme)
        self.enzyme = enzyme
        self.pattern = pattern
        self.is_regex = is_regex
        self.with_revcomp = with_revcomp
        if self.with_revcomp and self.pattern:
            self.rev_pattern = reverse_complement(pattern)

    def __call__(self, sequence):

        if self.biopython_enzyme is not None:
            return self.biopython_enzyme.search(sequence) == []

        if self.is_regex:
            cm_pattern = re.compile(self.pattern)
            if cm_pattern.search(sequence) is not None:
                if self.with_revcomp:
                    sequence_rev = reverse_complement(sequence)
                    return cm_pattern.search(sequence_rev) is not None
                else:
                    return True
            else:
                return False
        else:
            if self.pattern not in sequence:
                if self.with_revcomp:
                    return self.rev_pattern not in sequence
                else:
                    return True
            else:
                return False

    def __repr__(self):
        return "No pattern '%s'" % (self.pattern)


class SequenceLengthConstraint:
    def __init__(self, min_length=0, max_length=None):
        self.min_length = min_length
        self.max_length = max_length

    def __call__(self, sequence):
        L = len(sequence)
        upper_bound = self.max_length if self.max_length is not None else L
        return self.min_length <= L <= upper_bound

    def __str__(self):
        left_side = "" if (self.min_length == 0) else ("%d < " % self.min_length)
        right_side = "" if (self.max_length is None) else (" < %d" % self.max_length)
        return left_side + "length" + right_side


class GcContentConstraint:
    def __init__(self, min_gc=0, max_gc=1.0, memoize=False):
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.memoize = True
        self.memoization_dict = {}

    def __call__(self, sequence):
        if self.memoize:
            if sequence not in self.memoization_dict:
                result = self.min_gc <= gc_content(sequence) <= self.max_gc
                self.memoization_dict[sequence] = result
            return self.memoization_dict[sequence]
        return self.min_gc <= gc_content(sequence) <= self.max_gc

    def __str__(self):
        left_side = (
            "" if (self.min_gc == 0) else ("%.01f" % (self.min_gc * 100) + "% < ")
        )
        right_side = (
            "" if (self.max_gc == 1) else (" < %.01f" % (self.max_gc * 100) + "%")
        )
        return left_side + "GC" + right_side
