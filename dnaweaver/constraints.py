from .biotools import gc_content, reverse_complement
import re

class NoPatternConstraint:
    """Class of callables (sequence)-> True/False whether the sequence contains
    the pattern.

    Can be useful for defining constraints in DNA assembly methods or
    DNA providers.

    The interest of having this as a class is that a DnaSource using this
    constraint can be displayed as a string with the pattern appearing
    explicitly, which would not be the case for a function
    """

    def __init__(self, pattern, is_regex=False, with_revcomp=True):
        self.pattern = pattern
        self.is_regex = is_regex
        self.with_revcomp = with_revcomp

    def __call__(self, sequence):

        if self.is_regex:
            cm_pattern = re.compile(self.pattern)
            if (cm_pattern.search(sequence) is not None):
                if self.with_revcomp:
                    sequence_rev = reverse_complement(sequence)
                    return (cm_pattern.search(sequence_rev) is not None)
                else:
                    return True
            else:
                return False
        else:
            pattern_rev = reverse_complement(self.pattern)
            if self.pattern not in sequence:
                if self.with_revcomp:
                    return (pattern_rev not in sequence)
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
        left_side = ("" if (self.min_length == 0) else
                     ("%d < " % self.min_length))
        right_side =("" if (self.max_length is None) else
                     (" < %d" % self.max_length))
        return left_side + "length" + right_side


class GcContentConstraint:

    def __init__(self, min_gc=0, max_gc=1.0):
        self.min_gc = min_gc
        self.max_gc = max_gc

    def __call__(self, sequence):
        return self.min_gc <= gc_content(sequence) <= self.max_gc

    def __str__(self):
        left_side = ("" if (self.min_gc == 0) else
                     ("%.01f" % (self.min_gc*100) + "% < "))
        right_side =("" if (self.max_gc == 1) else
                     (" < %.01f" % (self.max_gc*100) + "%"))
        return left_side + "GC" + right_side

class PerBasepairPricing:

    def __init__(self, per_basepair_price):
        self.per_basepair_price = per_basepair_price
        self.min_basepair_price = per_basepair_price

    def __call__(self, sequence):
        return len(sequence) * self.per_basepair_price

    def __str__(self):
        return "$%.03f/bp" % self.per_basepair_price
