# Attempt import of optional module DNA Chisel
try:
    from dnachisel import Specification, DnaOptimizationProblem

    DNACHISEL_AVAILABLE = True
except:
    Specification = type("BLANK")  # fake type that matches nothing
    DNACHISEL_AVAILABLE = False


class ConstraintsMixin:
    def verify_constraints(self, sequence):
        """Return True iff `sequence` passes all `self.sequence_constraints`.

        Will automatically process DNA-Chisel constraints that would be in
        `self.sequence_constraints`.
        """
        constraints = self.sequence_constraints
        if not hasattr(self, "dnachisel_constraints"):
            self.dnachisel_constraints = [
                constraint
                for constraint in self.sequence_constraints
                if isinstance(constraint, Specification)
            ]

        if self.dnachisel_constraints != []:
            if not DNACHISEL_AVAILABLE:
                raise ImportError(
                    "Spotted DNA Chisel constraints, while "
                    "DNA Chisel is not installed."
                )
            # We provide an empty mutation space so it won't be recomputed
            # (which would take time and is useless here!)
            problem = DnaOptimizationProblem(
                sequence, self.dnachisel_constraints, mutation_space=[]
            )
            constraints = [
                constraint
                for constraint in constraints
                if not isinstance(constraint, Specification)
            ] + [lambda seq: problem.all_constraints_pass()]

        return all(constraint(sequence) for constraint in constraints)
