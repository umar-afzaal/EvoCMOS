# run.py (github.com/umar-afzaal/)
# This file is part of EvoMos
# This file contains the main routine


from optimizer import Optimizer


if __name__ == "__main__":
    
    # Configuration
    target = 'compressor_34T.v'
    graphCols = 46
    graphRows = 1
    printAfter = 100
    kmax = 2000000
    popSize = 20
    pa = 0.3
    bitsToMutate = 1
    seedPop = True

    # Configure and run the optimizer
    optimizer = Optimizer(kmax,
                          pa,
                          bitsToMutate,
                          graphCols,
                          graphRows,
                          popSize,
                          target,
                          printAfter,
                          seedPop)
    optimizer.run()

