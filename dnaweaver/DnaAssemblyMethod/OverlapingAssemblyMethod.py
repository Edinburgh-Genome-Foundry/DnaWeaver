from .DnaAssemblyMethod import DnaAssemblyMethod


class OverlapingAssemblyMethod(DnaAssemblyMethod):
    """General class for all overlaping assembly methods.

    Parameters
    ----------

    homology_arm_length
      Length of the homology arm, or "overhang". A length of L means that
      consecutive segments will overlap by 2*L

    """

    name = "Overlaping Assembly"

    def __init__(self, overhang_selector, **properties):
        super(OverlapingAssemblyMethod, self).__init__(**properties)
        selector = overhang_selector
        self.overhang_selector = selector
        if selector.has_location_filter:
            self.cut_location_constraints.append(
                selector.location_filter_method
            )
        self.compute_fragment_for_sequence_segment = (
            selector.compute_fragment_for_sequence_segment
        )

    # def extend_sequence(self, sequence):
    #     """Extend the end sequence with a start homology for circular assembly.
    #     """
    #     if self.topology == "circular":
    #         selector = self.overhang_selector
    #         start, end = selector.compute_segment_location(sequence, 0)
    #         size = end - start
    #         if size % 2:
    #             size += 1
    #         half_size = size // 2
    #         return sequence[-half_size:] + sequence + sequence[:half_size]
    #     else:
    #         return sequence


class GibsonAssemblyMethod(OverlapingAssemblyMethod):
    """Gibson Assembly Method. Just another overlap-method"""

    name = "Gibson Assembly"


class BuildAGenomeAssemblyMethod(OverlapingAssemblyMethod):
    """The Build-a-Genome Assembly Method. Just another overlap-method"""

    name = "Build-a-Genome"
