
class AssemblyMethod:
    pass


class GibsonAssemblyMethod(AssemblyMethod):

    def __init__(self, overlap=40):
        self.overlap = overlap

    def compute_fragments_sequences(self, cuts, sequence):
        L = len(sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        zones = zip(cuts, cuts[1:])
        overlap = self.overlap
        return {
            zone: sequence[max(0, zone[0] - overlap):
                           min(L, zone[1] + overlap)]
            for zone in zones
        }

    def validate_cuts(self, cuts, sequence):
        return True
