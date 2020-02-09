from .DnaAssemblyMethod import DnaAssemblyMethod


class BluntEndAssemblyMethod(DnaAssemblyMethod):
    def compute_fragment_for_sequence_segment(self, sequence, segment, **kw):
        start, end = segment
        return sequence[start:end]
