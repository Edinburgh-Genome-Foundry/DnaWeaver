class DnaAssemblyMethod(object):
    """General class for assembly methods.

    All assembly methods can store the following attributes:

    Parameters
    ----------

    duration
      Duration required for the assembly (e.g. an estimated upper bound). The
      time unit is left at the choice of the user.

    cut_location_constraints
      List or tuple of functions `(sequence, int) -> bool` which return for
      a sequence and an index (cut location) whether the location should be
      considered as a cutting site compatible with this assembly method (True)
      or not.
      The locations considered in fine are the locations which pass every
      constraint in the `cut_location_contraints` list.

    segment_constraints
      List or tuple of functions `(sequence, start, end) -> bool` which returns
      for a subsegment `(start, end)` whether the segment is a valid segment
      for the assembly.
      The segments considered in fine are the segments which pass every
      constraint in the `segments_constraints` list.

    min_segment_length
      Minimal length of the fragments that this assembly method allows

    max_segment_length
      Maximal length of the fragments that this assembly method allows


    force_cuts
      A function `sequence->[i1, i2, i3...]` which for a given sequence returns
      forced cuts locations.

    max_fragments
       Maximal number of fragments which can be assembled with this method.

    sequence_constraints
      General constraints on the sequence so that it can be built with this
      assembly method. It is a list of functions `seq->bool`. If one of
      these constraints returns False the sequence is refused.

    reference
      Reference to e.g. a paper or a protocol describing the method.
    """

    name = "None"

    def __init__(
        self,
        duration=0,
        cost=0,
        reference=None,
        cut_location_constraints=(),
        segment_constraints=(),
        min_segment_length=0,
        max_segment_length=None,
        force_cuts=(),
        suggest_cuts=(),
        max_fragments=None,
        sequence_constraints=(),
        cuts_set_constraints=(),
    ):
        self.duration = duration
        self.cost = cost
        self.cut_location_constraints = list(cut_location_constraints)
        self.segment_constraints = list(segment_constraints)
        self.min_segment_length = min_segment_length
        self.max_segment_length = max_segment_length
        self.sequence_constraints = list(sequence_constraints)
        self.max_fragments = max_fragments
        self.cuts_set_constraints = list(cuts_set_constraints)
        self.reference = reference

        if callable(suggest_cuts):
            self.suggest_cuts = suggest_cuts
        else:
            # means the cuts are a constant
            self.suggest_cuts = lambda *a: suggest_cuts

        if callable(force_cuts):
            self.force_cuts = force_cuts
        else:
            # means the forced cuts are a constant
            self.force_cuts = lambda *a: force_cuts

    def dict_description(self):
        result = {
            "name": self.name,
            "minimum segment length": self.min_segment_length,
            "maximum segment length": self.max_segment_length,
            "maximal number of fragments": self.max_fragments,
            "reference": self.reference,
            "duration": self.duration,
            "cost": self.cost,
        }
        result.update(self.additional_dict_description())
        return result

    # def extend_sequence(self, sequence):
    #     return sequence

    def additional_dict_description(self):
        return {}
