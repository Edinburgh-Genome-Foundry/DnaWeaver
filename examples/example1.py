from dnaadvisor import optimize_cuts_with_GA, DnaOffer, DnaOrderingProblem, GibsonAssemblyMethod
from dnachisel import random_dna_sequence, PROVIDERS_CONSTRAINTS
sequence = random_dna_sequence(200) + 'GGTCTC' + random_dna_sequence(300) + 'GGTCTC' + random_dna_sequence(200)
sequence = 'GGTCTC' +  random_dna_sequence(200)
gen9_offer = DnaOffer(
    name="Gen9Offer",
    constraints = PROVIDERS_CONSTRAINTS["Gen9"],
    pricing = lambda sequence: 0.08*len(sequence)
)

idt_offer = DnaOffer(
    name="IDTOffer",
    constraints= PROVIDERS_CONSTRAINTS["IDT"],
    pricing = lambda sequence: 0.10*len(sequence)
)

twix_offer=DnaOffer(
    name="Twix",
    constraints = [],
    pricing = lambda sequence: 0.20*len(sequence)

)

prob = DnaOrderingProblem(
    sequence= sequence,
    offers = [gen9_offer, idt_offer, twix_offer],
    assembly_method= GibsonAssemblyMethod(20)
)
pop, log, stats, hof, history = optimize_cuts_with_GA(prob, 5, pop_size=80, generations=300)
winner = list(hof[0])
print prob.score_cuts(winner)
offers = prob.find_best_offers(winner).values()
offers = sorted(offers, key=lambda o: o.zone)
for offer in offers:
    print offer

def random_constrained_sequence(length, constraints):
    dna = random_dna_sequence(300)
    canvas = DnaCanvas(sequence, constraints=constraints)
    canvas.solve_all_constraints_one_by_one()
    return canvas.sequence
