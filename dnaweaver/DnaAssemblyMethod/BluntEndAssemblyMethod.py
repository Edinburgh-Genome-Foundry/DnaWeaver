from .DnaAssemblyMethod import DnaAssemblyMethod


class BluntEndAssemblyMethod(DnaAssemblyMethod):
    def compute_sequence_fragment(self, sequence, segment):
        start, end = segment
        return sequence[start:end]
