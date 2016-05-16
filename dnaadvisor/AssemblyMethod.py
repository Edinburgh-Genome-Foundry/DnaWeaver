from dnachisel.biotools import reverse_complement

class AssemblyMethod:
    """General class for assembly methods.

    Yeah that class is useless right now but bear with me.
    """
    pass


class OverlapingAssemblyMethod(AssemblyMethod):
    """Gibson Assembly Method.

    Parameters
    ----------

    homology_arm_length
      Length of the homology arm, or "overhang". A length of L means that
      consecutive segments will overlap by 2*L

    """

    def __init__(self, homology_arm_length=20):
        self.homology_arm_length = homology_arm_length


    def compute_fragment_sequence(self, segment, sequence):
        """Return the segment's sequence with flanking sequences.

        Parameters
        ----------

        segment
          A pair of integers (start, end) delimiting the subfragment
          sequence[start:stop]

        sequence
          An "ATGC" DNA sequence string

        """
        L = len(sequence)
        start, end = segment
        return sequence[max(0, start - self.homology_arm_length):
                        min(L, end + self.homology_arm_length)]

    def compute_fragments_sequences(self, cuts, sequence):
        """Compute the sequences (with flanks) of all fragments corresponding
        to the cuts.

        Parameters
        ----------

        cuts
          a list of cuts between 0 and len(sequence) (if these two numbers are
          omited they will be added anyways)

        sequence
          An "ATGC" DNA sequence string

        """
        L = len(sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        return {
            segment: self.compute_fragment_sequence(segment, sequence)
            for segment in zip(cuts, cuts[1:])
        }

class GibsonAssemblyMethod(OverlapingAssemblyMethod):
    """Gibson Assembly Method. Just another overlap-method"""

class BuildAGenomeAssemblyMethod(OverlapingAssemblyMethod):
    """The Build-a-Genome Assembly Method. Just another overlap-method"""

class GoldenGateAssemblyMethod(AssemblyMethod):
    """The Golden Gate Assembly Method.

    This method adds overhangs with a Type IIS REase site to the segments.

    Parameters
    ----------

    left_overhang
      An ATGC DNA sequence representing the left overhang, in 5' -3'.
      For practicality, this can contain "[BsaI]", "[BsmBI]", "[BbsI]", which
      will be replaced by their restriction sites


    right_overhang
      An ATGC DNA sequence representing the right overhang, in 5'-3'.
      Be careful, I said 5'-3' !!!
      If left to None, will be equal to left_overhang.

    Examples
    --------

    >>> #  The fragments will have an overhang with an homology region, a BsaI
    >>> #  site, and a "T" for the wildcard
    >>> assembly_method = GoldenGateAssemblyMethod(
    >>>  left_overhang = "ATGTGTCGTGTGTGCGTA[BsaI]T",
    >>>  right_overhang = "TTCTCTCGATAAATGGCC[BsaI]T"
    >>> )

    """
    enzymes_dict = {
        "[BsaI]": "GGTCTC",
        "[BsmBI]": "CGTCTC",
        "[BbsI]": "GAAGAC",
    }

    def __init__(self, left_overhang="[BsaI]A", right_overhang=None):
        if right_overhang is None:
            right_overhang = left_overhang
        for name, site in self.enzymes_dict:
            left_overhang = left_overhang.replace(name, site)
            right_overhang = right_overhang.replace(name, site)

        self.left_overhang = left_overhang
        self.right_overhang = right_overhang
        self.right_overhang_rev = reverse_complement(right_overhang)


    def compute_fragment_sequence(self, segment, sequence):
        """Return the segment's sequence with flanking sequences for

        Parameters
        ----------

        segment
          A pair of integers (start, end) delimiting the subfragment
          sequence[start:stop]

        sequence
          An "ATGC" DNA sequence string

        """
        L = len(sequence)
        start, end = segment
        segment = sequence[max(0, start - 2): min(L, end + 2)]
        return self.left_overhang + segment + self.right_overhang_rev

    def compute_fragments_sequences(self, cuts, sequence):
        """Compute the sequences (with flanks) of all fragments corresponding
        to the cuts.

        Parameters
        ----------

        cuts
          a list of cuts between 0 and len(sequence) (if these two numbers are
          omited they will be added anyways)

        sequence
          An "ATGC" DNA sequence string

        """
        L = len(sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        return {
            segment: self.compute_fragment_sequence(segment, sequence)
            for segment in zip(cuts, cuts[1:])
        }
