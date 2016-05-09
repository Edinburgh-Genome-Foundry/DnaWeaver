"""Optimization techniques"""

import numpy as np
from deap import creator, base, tools, algorithms
from tqdm import tqdm
import itertools as itt
import networkx as nx

def optimize_cuts_with_graph(sequence, segment_score_function,
                             segment_length_range=(500, 2000),
                             nucleotide_resolution=10):
    min_length, max_length = segment_length_range
    nodes = range(0, len(sequence), nucleotide_resolution)
    nodes[-1] = len(sequence)
    segments = [
        (start, end)
        for start, end in itt.product(nodes, nodes)
        if min_length < end - start < max_length
    ]
    graph = nx.DiGraph()
    for start, end in tqdm(segments):
        cost = segment_score_function((start, end), sequence)
        graph.add_edge(start, end, weight=cost)
    return nx.dijkstra_path(graph, 0, len(sequence))


def optimize_costs_with_graph(dna_ordering_problem, length_range=(500, 2000),
                              nucleotide_resolution=10):
    def segment_score_function(segment, sequence):
        offers = dna_ordering_problem.find_offers(segment)
        if offers == []:
            return float(1e8)
        else:
            return min([offer.price for offer in offers])
    cuts = optimize_cuts_with_graph(
        dna_ordering_problem.sequence,
        segment_score_function,
        length_range=length_range,
        nucleotide_resolution=nucleotide_resolution
    )
    return dna_ordering_problem.find_best_offers(cuts)


def optimize_cost_with_GA(dna_ordering_problem, ncuts, pop_size=20,
                          generations=40):
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))

    creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    L = len(dna_ordering_problem.sequence)
    toolbox.register("attr_int", np.random.randint, 0, L)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attr_int, ncuts)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    def evaluate(individual):
        return -dna_ordering_problem.score_cuts(list(set(individual))),

    mutation_amplitude = L / ncuts

    def mutate(individual):
        individual += np.random.randint(-mutation_amplitude,
                                        mutation_amplitude, ncuts)
        individual[individual < 0] = 0
        individual[individual > L] = L

        return individual,

    def mutate(individual):
        fitness = individual.fitness
        individual += np.random.randint(-mutation_amplitude,
                                        mutation_amplitude, len(individual))
        individual = np.minimum(L, np.maximum(0, individual))
        individual = creator.Individual(np.array(list(set(individual))))
        individual.fitness = fitness
        return individual,

    toolbox.register("evaluate", evaluate)
    toolbox.register("mutate", mutate)
    toolbox.register("select", tools.selTournament, tournsize=3)

    history = tools.History()
    toolbox.decorate("mutate", history.decorator)

    hof = tools.HallOfFame(1, similar=np.array_equal)
    stats = tools.Statistics()
    #stats.register("mean", lambda ind: np.mean(ind))
    pop = toolbox.population(n=pop_size)
    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0, mutpb=0.2,
                                   ngen=generations, stats=stats,
                                   halloffame=hof, verbose=False)
    return pop, log, stats, hof, history
