
class AssemblyMethod:
    pass


class GibsonAssemblyMethod(AssemblyMethod):

    def __init__(self, overlap=40):
        self.overlap = overlap


    def compute_fragment_sequence(self, zone, sequence):
        L = len(sequence)
        return sequence[max(0, zone[0] - self.overlap):
                        min(L, zone[1] + self.overlap)]

    def compute_fragments_sequences(self, cuts, sequence):
        L = len(sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        return {
            zone: self.compute_fragment_sequence(zone, sequence)
            for zone in zip(cuts, cuts[1:])
        }
