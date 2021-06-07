#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H

#include <random>

// Method of determing selection of survivors that will carry properties into the new individuals of the newGeneration
// Called from optimize::optimization.cu
// Input: pool - (inputParameters) generation of individuals to pull from, no implied order
//        poolSize - integer number of length of the pool
//        selectionSize - integer number of how many survivors to choose out of the pool
//        survivors - pointer array of individuals to promote via crossover / mutate 
// Output: pool is left sorted by individuals with higher velocity difference
//         survivors contains an array of size selectionSize of individuals to be used in newGeneration()
void selectSurvivors(Individual * pool, int poolSize, int selectionSize, Individual* survivors, const double & ratio);

///////////////////////////////////////////////////////////////
// Crossover Functions                                       //
///////////////////////////////////////////////////////////////

//Enumeration Documentation
// PARTNER1 (1) - At this index of the OPTIM_VARS array, accept the gene value from parent 1
// PARTNER2 (2) - At this index of the OPTIM_VARS array, accept the gene value from parent 2
//      AVG (3) - At this index of the OPTIM_VARS array, use the average value for the gene from parent 1 and parent 2

// Mask definition / example
// mask is a pointer array set to OPTIM_VARS length and contains enumeration values (ranges from 1 to 3)
// used in determing for a given parameter value which parent to receive it from
// Example: (with OPTIM_VARS being 6)
//          [1,1,1,1,1,1] - the individual accepts all gene values from parent 1, effectively being a copy
//          [1,1,1,2,2,2] - the first half of the parameter array is to be copied from parent 1, the other half is to be copied from parent 2
//          [3,3,3,3,3,3] - for all gene values, the new individual contains the average between the two parents gene values
// flipping the mask is changing the 1 values to 2, and 2 values to 1 to create a mask that is opposite to maximize number of new individuals created with crossover
// average is not changed when flipping the mask as there is no opposite equivalent (instead relying on mutation for the new individuals to be different)

// Creates a random bifurcation mask
// ** currently not in use - replaced with bundleVars / average **
// Randomly picks one index to be the start of the '2's from mask
void crossOver_randHalf(int * mask, std::mt19937_64 & rng);


// Creates a mask where no mixing occurs
// ** currently not in use **
// Sets the mask to be entirely PARTNER1, use flipMask() to have PARTNER2
void crossOver_oneParent(int * mask);

// Creates a mask that is contains randomly chosen values in each index
// Each element in a mask is randomly set to either PARTNER1 or PARTNER2
void crossOver_wholeRandom(int * mask, std::mt19937_64 & rng);

// Generates crossover mask that maintains paramter relationships (gamma, tau, coast values grouped)
// Similar to crossOver_wholeRandom, but with parameter grouping
// Input: mask - set, size OPTIM_VARS, 
//        rng - managed from optimize()
// Output: mask contains values either 1 or 2 (equivalent to PARTNER1 or PARTNER2)
void crossOver_bundleVars(int * mask, std::mt19937_64 & rng);

// Sets the entire mask to be AVG for length OPTIM_VARS
// Input: mask - pointer integer array of length OPTIM_VARS
// Output: mask is set to contain all AVG values
void crossOver_average(int * mask);

// Utility to flip the polarity of a mask
// Input:  mask should have already be set from other crossover functions
// Output: each PARTNER1 in mask will be reassigned to be a PARTNER2
//         each PARTNER2 will be reassigned PARTNER1
//         AVG is left unchanged
void flipMask(int * mask);

// Copy contents of maskIn into maskOut of size OPTIM_VARS
// Output: maskIn remains unchanged, maskOut will have same contents as maskIn
void copyMask(int * maskIn, int * maskOut);

///////////////////////////////////////////////////////////////////////////////

// Utility function for mutate() to get a random double with high resolution
// Input: max - the absolute value of the min and max(min = -max) of the range
// Output: A double value that is between -max and +max
double getRand(double max, std::mt19937_64 & rng);

// Creates a new rkParameters individual by combining properties of two parent Individuals using a crossover mask
// Input: two rkParameter individuals (p1 and p2) - source of genes for new individual
//        mask - Contains maskValue values and length of OPTIM_VARS,
//               determines how the two parent properties are merged into creating the new individual
//        cConstants, annealing, rng, generation - passed through to mutate()
// Output: Returns rkParameter object that is new individual
// Called from generateChildrenPair, calls mutate
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int * mask, const cudaConstants * cConstants, double annealing, std::mt19937_64 & rng, double generation);

// Utility function, generates a boolean mask for which paramters to mutate (1: mutate, 0: not mutated)
// Number of genes mutated is a compound probability of n-1 genes before it
// First gene chance is base mutation rate, second chance is after first is mutated
// chance of mutating n genes = mutation_rate^n
// input: rng - random number generating object, defined in optimize()
//        mutateMask - pointer to a boolean array, length of OPTIM_VARS
//                   - sets genes being mutated to true and others to false
//        mutation_rate - a double value less than 1 that is the chance a gene will be mutated
//                      - called iteratively to mutate more genes
// output: mutateMask contains false for genes that are not mutating, true for genes that are to be mutated
// Called by mutate()
void mutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate);

// Handles potential mutation of individual 
// Calls mutateMask, then applied mutations as necessary
// mutate a gene by adding or subtracting a small, random value on a parameter property
// Input: p1 - rkParameter that is taken to be the mutation base
//        rng - random number generator to use
//        annealing - a scalar value on the max random number when mutating
//        cConstants - holds properties to use such as mutation rates and mutation scales for specific parameter property types
// Output: Returns rkParameter object that is the mutated version of p1
// Called by generateNewIndividual
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, const cudaConstants* gConstant, double generation);

// Method that creates a pair of new Individuals from a pair of parent individuals and a mask
// Input: pool - (output) pointer array to Individuals that is where the new pair of individuals are stored
//        survivors - (parents) pointer array to Individuals to access the two parents from
//        mask - pointer array of maskValues used to decide on which property from which parent is acquired (or average of the two)
//        newIndCount - value that tracks number of newly created indiviudals in the pool so far in the newGeneration process
//                    - also impacts where to put the new individuals in the pool
//        parentsIndex - value determing where the pair of parent survivors are selected (parent 1 is at parentsIndex, parent 2 is offset by +1)
//        annealing - double variable passed onto mutateNewIndividual
//        poolSize - length of the pool array
//        rng - random number generator passed on to generateNewIndividual
//        cConstants - passed on to generateNewIndividual
// Output: pool contains two newly created individuals at (poolSize - 1 - newIndCount) and (poolSize - 2 - newIndCount)
//         mask is flipped in polarity between each (refer to flipMask method) 
//         newIndCount is incremented by +2
// Called by newGeneration()
void generateChildrenPair(Individual *pool, Individual *survivors, int * mask, int& newIndCount, int parentsIndex, double annealing, int poolSize, std::mt19937_64 & rng, const cudaConstants* cConstants, double generation);

// Creates the next pool to be used in the optimize function in opimization.cu
// Input: survivors - (parents) Individual pointer array of Individuals to be used in creating new individuals
//        pool - Sorted Individual pointer array that contains current generation and is modified for next generation
//        survivorSize - length of survivors array
//        poolSize - length of pool array
//        annealing - passed onto generateChildrenPair
//        rng - Random number generator to use for making random number values
//        cConstants - passed onto generateChildrenPair
// Output: lower (survivorSize * 4) portion of pool is replaced with new individuals
//         Each parent pair produces 8 new children (4 masks, two children per mask)
//         Returns number of new individuals created (newIndCount)
int newGeneration(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, const cudaConstants* cConstants, std::mt19937_64 & rng, double generation );

#include "ga_crossover.cpp"
#endif
