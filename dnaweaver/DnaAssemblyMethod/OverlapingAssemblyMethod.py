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
        self.compute_sequence_fragment = selector.compute_sequence_fragment


class GibsonAssemblyMethod(OverlapingAssemblyMethod):
    """Gibson Assembly Method. Just another overlap-method"""

    name = "Gibson Assembly"


class BuildAGenomeAssemblyMethod(OverlapingAssemblyMethod):
    """The Build-a-Genome Assembly Method. Just another overlap-method"""

    name = "Build-a-Genome"
