# optimizer.py (github.com/umar-afzaal/)
# This file is part of EvoMos
# This file contains class definition for Optimizer


from population import Population
from chromosome import Chromosome
import os
import math
from operator import itemgetter, attrgetter
from random import choices, randint, uniform, choice           
from copy import deepcopy
import signal
from datetime import datetime
#from shutil import rmtree


class Optimizer():
    def __init__(self, kmax, pa, bitsToMutate, graphCols, graphRows, popSize, target, printAfter=100, seedPop=False, printBest=True):
        '''
        - param: <int> kmax: Maximum number of generations to run (Termination criteria)
        - param: <float> pa: Parasitic probability
        '''
        self.pop = Population(graphCols, graphRows, popSize, target, seedPop)
        self.kmax = kmax
        self.nCuckooEggs = int(popSize*pa)
        self.kway = int(popSize*pa)
        self.printAfter = printAfter
        self.printBest = printBest
        self.h = bitsToMutate

        # print optimizer configuration to console
        currentTime = datetime.now()
        datetimeString = currentTime.strftime("%d/%m/%Y %H:%M:%S")
        print(f"Timestamp: {datetimeString}")
        print(f"Target: {target}")
        print(f"Graph: {graphCols}x{graphRows}")
        print(f"Pop size: {popSize}")
        print(f"pa: {pa}")
        print(f"kmax: {kmax}")
        print(f"Step size (bits): {bitsToMutate}")
        print(f'Print after: {printAfter}\n\n')


    def sig_int_handler(self, signal, frame):  # Press Ctrl+C to interrupt
        print ("\n\nABORTED!")
        try:
            os.remove("dsn")
            os.remove("tb_ckt.v")
            os.remove("ckt.v")
            os.remove("ckt.vcd")
            #os.rmtree("__pycache__")
        except: pass
        try:
            self.pop.sort()
            print(f"#{self.k}",
                  f"(f:{round(self.pop.population[0].fitness, 3)},",
                  f"er:{round(self.pop.population[0].er, 3)},",
                  f"hd:{round(self.pop.population[0].hd, 3)})",
                  f"(cells:{self.pop.population[0].gates})",
                  f"(ckts:{len(self.front)})")
            if self.printBest == True:
                print ("__________________")
                for line in self.pop.population[0].netlist:
                    print(line.strip())
            self.dump_to_file()
        except: pass
        os._exit(0)


    def update(self):
        if self.k%self.printAfter==0:
            print(f"#{self.k}",
                  f"(f:{round(self.pop.population[0].fitness, 3)},",
                  f"er:{round(self.pop.population[0].er, 3)},",
                  f"hd:{round(self.pop.population[0].hd, 3)})",
                  f"(cells:{self.pop.population[0].gates})",
                  f"(ckts:{len(self.front)})")
            if self.printBest==True:
                print ("__________________")
                for line in self.pop.population[0].netlist:
                    print(line.strip())
                print('')
                
        if self.printSwitch == True:
            if self.pop.population[0].er==0.0:
                print(f'Functional-correctness achieved at #{self.k}')
                self.printSwitch = False
                print(f"#{self.k}",
                      f"(f:{round(self.pop.population[0].fitness, 3)},",
                      f"er:{round(self.pop.population[0].er, 3)},",
                      f"hd:{round(self.pop.population[0].hd, 3)})",
                      f"(cells:{self.pop.population[0].gates})",
                      f"(ckts:{len(self.front)})")
                if self.printBest==True:
                    print ("__________________")
                    for line in self.pop.population[0].netlist:
                        print(line.strip())
                    print('')


    def dump_to_file(self):
        with open('csv', 'w') as f:
            f.write('gates*netlist\n')
            for nmos,pmos,gates in sorted(self.front.keys(), key=itemgetter(2)): # itemgetter(2) sorts by gates
                f.write( f'{gates}'+'*'+(''.join(map(str, self.front[(nmos,pmos,gates)])).replace('\n','\\n'))+'\n' )


    def offsprings(self):
        father = self.pop.population[0] # Best member as father for breeding
        children = []
        
        while(True):
            # Tournament Selection for mother
            tournament = choices(self.pop.population[1:], k=self.kway) # eggs for tournament
            mother = min(tournament, key=attrgetter("fitness")) # winner

            ## Single-point Crossover
            iBit = randint(0,father.Lg-1)

            child1 = Chromosome()
            child1.set_params(father.ni, father.no, father.nc, father.nr)
            
            child2 = Chromosome()
            child2.set_params(father.ni, father.no, father.nc, father.nr)

            if uniform(0,1) < uniform(0,1):
                child1.chromosome = father.chromosome[:iBit]+mother.chromosome[iBit:]
                child2.chromosome = mother.chromosome[:iBit]+father.chromosome[iBit:]
            else:
                child1.chromosome = mother.chromosome[:iBit]+father.chromosome[iBit:]
                child2.chromosome = father.chromosome[:iBit]+mother.chromosome[iBit:]

            child1.point_mutation(child1.chromosome, h=self.h) # Mutation
            child1.calc_fitness(self.pop.seed.simulation)
            
            child2.point_mutation(child2.chromosome, h=self.h) # Mutation
            child2.calc_fitness(self.pop.seed.simulation)

            children.append(child1)
            children.append(child2)

            if len(children) >= self.nCuckooEggs:
                return( children[:self.nCuckooEggs] )

            
    def run(self):
        signal.signal(signal.SIGINT, self.sig_int_handler)

        # create levy flight steps array
##        self.levy = []
##        for i in range(1000):    
##            self.levy.append( round(math.pow(uniform(0.0001,0.9999),-1.0/3.0)) )

        # search loop
        counter=0
        self.printSwitch = True
        self.front = {}
        for k in range(self.kmax+1):
            # lay cuckoo eggs
            for i in range(self.nCuckooEggs):
                index = randint(0,self.pop.popSize-1)
                cuckooEgg = deepcopy(self.pop.population[index]) # Pick a member of population at random
                cuckooEgg.point_mutation(cuckooEgg.chromosome, h=self.h)
                cuckooEgg.calc_fitness(self.pop.seed.simulation)

                if (cuckooEgg.fitness <= self.pop.population[index].fitness): # Check if the cuckoo egg is better
                    self.pop.population[index] = cuckooEgg

            self.pop.sort() # sort population

            if counter<25:
                # genetically modify eggs with probability pa
                if self.nCuckooEggs > 0:
                    children = self.offsprings()
                    for i in range(self.nCuckooEggs):
                        self.pop.population[(self.pop.popSize-1)-(i)] = children[i]
                counter+=1
            else:
                counter=0
                for i in range(self.nCuckooEggs):
                    self.pop.population[(self.pop.popSize-1)-(i)] = self.pop.get_random_candidate()

            # store functionally-correct circuits
            funCorrect = [egg for egg in self.pop.population if egg.er==0.0] # Consider only functionally-correct
            if funCorrect != []: # If there are some
                for egg in funCorrect:
                    if (egg.nmos, egg.pmos, egg.gates) not in self.front.keys(): # Check if already stored
                        self.front[(egg.nmos, egg.pmos, egg.gates)] = egg.netlist

            # status update
            self.k=k
            self.update()

        # dump collected circuits to file
        self.dump_to_file()



        
