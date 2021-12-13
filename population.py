# population.py ( github.com/umar-afzaal/ )
# This file is part of EvoMos
# This file contains definition of the Population class


from chromosome import Chromosome
from operator import attrgetter


class Population():
    def __init__(self, nc, nr, popSize, target, seedPop=False):
        # Create seed object of class Chromosome
        seed = Chromosome()
        if '.pla' in target:
            seedPop=False
            seed.read_truth_table(target)
        else:        
            if seedPop==True:
                popSize = popSize-1
                seed.encode(nc, nr, target)
                seed.decode()
            else:
                seed.read_netlist(target)
                seed.set_params(seed.ni, seed.no, nc, nr)
            seed.simulate()
            seed.calc_fitness(seed.simulation)

        # Create the population    
        self.popSize = popSize
        self.population = []

        # Rest of the population
        for i in range(popSize):
            egg = Chromosome()
            egg.set_params(seed.ni, seed.no, nc, nr)
            egg.get_random_chromosome()
            egg.calc_fitness(seed.simulation)
            self.population.append(egg)
            
        # Append seed to population
        if seedPop==True:
            self.population.append(seed)
            
        self.seed = seed
        self.nc = nc
        self.nr = nr


    def sort(self):
        # sort the population in ascending order of error
        self.population = sorted(self.population, key=attrgetter("fitness"))


    def get_random_candidate(self):
        egg = Chromosome()
        egg.set_params(self.seed.ni, self.seed.no, self.nc, self.nr)
        egg.get_random_chromosome()
        egg.calc_fitness(self.seed.simulation)
        return (egg)


