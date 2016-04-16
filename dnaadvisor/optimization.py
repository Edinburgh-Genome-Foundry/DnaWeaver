"""Optimization techniques"""

import numpy as np
from deap import creator, base, tools, algorithms


def optimize_cuts_with_GA(dna_ordering_problem, ncuts, pop_size=20, generations=40):
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
