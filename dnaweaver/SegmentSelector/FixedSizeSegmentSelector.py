from .SegmentSelector import SegmentSelector


class FixedSizeSegmentSelector(SegmentSelector):
    """Selects segments of a constant size.

    Great for methods involving large homology regions where melting
    temperature matters less.
    """

    def __init__(self, segment_size=100, left_addition="", right_addition=""):
        self.segment_size = segment_size
        self.left_addition = left_addition
        self.right_addition = right_addition

    def compute_segment_location(self, sequence, index):
        return self.get_segment_coordinates(index, self.segment_size, len(sequence))

    @property
    def max_homology_size(self):
        return self.segment_size

    def __str__(self):
        result = "FixedSize(%dbp)" % self.segment_size
        if self.left_addition:
            result = ("...%s-" % self.left_addition[-12:]) + result
        if self.right_addition:
            result = result + ("-%s..." % self.right_addition[:12])
        return result
