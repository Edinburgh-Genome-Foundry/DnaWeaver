from .DnaAssemblyMethod import DnaAssemblyMethod
from ..biotools import reverse_complement


class OverlapingAssemblyMethod(DnaAssemblyMethod):
    """General class for all overlapping assembly methods.

    Parameters
    ----------

    homology_arm_length
      Length of the homology arm, or "overhang". A length of L means that
      consecutive segments will overlap by 2*L.
    """

    name = "Overlaping Assembly"
    alternate_fragments_orientation = False

    def __init__(self, overhang_selector, **properties):
        super(OverlapingAssemblyMethod, self).__init__(**properties)
        selector = overhang_selector
        self.overhang_selector = selector
        if selector.has_location_filter:
            self.cut_location_constraints.append(selector.location_filter_method)

    def compute_fragment_for_sequence_segment(self, sequence, segment, **kw):
        selector = self.overhang_selector.compute_fragment_for_sequence_segment
        fragment = selector(sequence, segment)
        if self.alternate_fragments_orientation:
            if kw.get("segment_position", 0) % 2:
                fragment = reverse_complement(fragment)
        return fragment


class GibsonAssemblyMethod(OverlapingAssemblyMethod):
    """Gibson Assembly Method. Just another overlap-method"""

    name = "Gibson Assembly"


class OligoAssemblyMethod(OverlapingAssemblyMethod):
    """The Build-a-Genome Assembly Method. Just another overlap-method"""

    alternate_fragments_orientation = True
    name = "Oligo Assembly"
