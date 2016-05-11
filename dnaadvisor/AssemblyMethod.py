
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
