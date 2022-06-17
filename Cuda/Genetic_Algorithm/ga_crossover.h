#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H

#include <random>
#include <chrono>
#include <vector>

#include "..\Config_Constants\constants.h"
#include "..\Genetic_Algorithm\adult.h"
#include "..\Genetic_Algorithm\child.h"

///////////////////////////////////////////////////////////////
// Crossover Functions                                       //
///////////////////////////////////////////////////////////////

//Enumeration Documentation
// PARTNER1 (1) - At this index of the OPTIM_VARS array, accept the gene value from parent 1
// PARTNER2 (2) - At this index of the OPTIM_VARS array, accept the gene value from parent 2
//      AVG (3) - At this index of the OPTIM_VARS array, use the average value for the gene from parent 1 and parent 2

//TODO: Document org: probably group / document everything together

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
void crossOver_randHalf(std::vector<int> & mask, std::mt19937_64 & rng);


// Creates a mask where no mixing occurs
// ** currently not in use **
// Sets the mask to be entirely PARTNER1, use flipMask() to have PARTNER2
void crossOver_oneParent(std::vector<int> & mask);

// Creates a mask that is contains randomly chosen values in each index
// Each element in a mask is randomly set to either PARTNER1 or PARTNER2
void crossOver_wholeRandom(std::vector<int> & mask, std::mt19937_64 & rng);

// Generates crossover mask that maintains paramter relationships (gamma, tau, coast values grouped)
// Similar to crossOver_wholeRandom, but with parameter grouping
// Input: mask - set, size OPTIM_VARS, 
//        rng - managed from optimize()
// Output: mask contains values either 1 or 2 (equivalent to PARTNER1 or PARTNER2)
void crossOver_bundleVars(std::vector<int> & mask, std::mt19937_64 & rng);

// Sets the entire mask to be AVG for length OPTIM_VARS
// Input: mask - pointer integer array of length OPTIM_VARS
// Output: mask is set to contain all AVG values
void crossOver_average(std::vector<int> & mask);

// Utility to flip the polarity of a mask
// Input:  mask should have already be set from other crossover functions
// Output: each PARTNER1 in mask will be reassigned to be a PARTNER2
//         each PARTNER2 will be reassigned PARTNER1
//         AVG is left unchanged
void flipMask(std::vector<int> & mask);

// Copy contents of maskIn into maskOut of size OPTIM_VARS
// Output: maskIn remains unchanged, maskOut will have same contents as maskIn
void copyMask(int * maskIn, int * maskOut);

///////////////////////////////////////////////////////////////////////////////

// Creates a new rkParameters individual by combining properties of two parent Individuals using a crossover mask
// Input: two rkParameter from individuals (p1 and p2) - source of genes for new individual
//        mask - Contains maskValue values and length of OPTIM_VARS,
//               determines how the two parent properties are merged into creating the new individual
//        cConstants, annealing, rng, generation - passed through to mutate()
// Output: Returns rkParameter object that is new individual
// Called from generateChildrenPair, calls mutate
rkParameters<double> generateNewChild(const rkParameters<double> & p1, const rkParameters<double> & p2, const std::vector<int> & mask, const cudaConstants * cConstants, const double & annealing, std::mt19937_64 & rng, const double & generation);

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
//        generation - passed in to report on mutations if desired
//        mutationScale - scalar for the mutation intensity, allows for control on the intensity of mutations between children and duplicates
// Output: Returns rkParameter object that is the mutated version of p1
// Called by generateNewIndividual
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, const double & annealing, const cudaConstants* cConstants, const double & generation, const double & mutationScale);

// Utility function for mutate() to get a random double with high resolution
// Input: max - the absolute value of the min and max(min = -max) of the range
// Output: A double value that is between -max and +max
double getRand(double max, std::mt19937_64 & rng);

// Method that creates a pair of new Children from a pair of parent Adults and a mask
// Input:  parent1 - the first parent that the children will draw parameters from
//         parent2 - the second parent that the children will draw parameters from
//         newChildren - a pointer to the newChildren array that is used to add the newly generated children to
//         mask - pointer vector of maskValues used to decide on which property from which parent is acquired (or average of the two)
//         annealing - double variable passed onto mutateNewIndividual
//         rng - random number generator passed on to generateNewChild
//         numNewChildren - Used within newGeneration to keep track of how many children have been generated; added to whenever a new child is generated
//         generation - passed on to GenerateNewChild... Will eventually be used in reporting mutations
//         cConstants - passed on to generateNewChildren... Used for calculating mutations and constructing children
// Output: The newChildren array will contain two newly generated children
//         mask is flipped in polarity (refer to flipMask method) 
//         numNewChildren is incremented by +2
// Called by newGeneration() each time a new mask is generated 

//TODO: Add const back to parent1 and parent2
void generateChildrenPair( Adult & parent1, Adult & parent2, Child * newChildren, const int & generateNum, std::vector<int> & mask, const double & annealing, std::mt19937_64 & rng, int & numNewChildren, const int & generation, const cudaConstants* cConstants);

// newGeneration will generate a mask and use random adults from oldAdults to generate new children
//      It will then simulate the newChildren.
//      Finally, it will convert the new children into adults and insert them into the newAdults vector
//      Note: the number of new children generated should ideally equal the number of threads on the gpu
//          this is set in the config file
// Inputs:  oldAdults - a vector of adults that will be used to pull random parents to generate children
//          newAdults - a vector of adults that will be filled by children that have been simulated and converted into adults
//          annealing - a double that will be passed on to other functions; it determines how much a child's parameters are mutated
//          generation - the generation that the algorithim is on; this will be passed on to other functions and is used for reporting
//          rng - a random number generator that is both passed on to other functions and used to pull random parents
//          cConstants - the config constants; it is passed on to other functions and is used for constructing children and mutating their parameters
// Outputs: newAdults will be filled with newly generated adults; it will be ready to be sorted and compared to other generations
// NOTE: This function is called at the beginning of each generation within optimize() 
void newGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants);

// This function will generate a certain number of children via the ga_crossover method
//      It will ensure that duplicate adults within the survivor pool do not create children
// Inputs:  parents - this adult vector holds the non-duplicate parent list
//          newChildren - this is the array that generated children will be inserted into
//          numNeededChildren - tracks how many children needs to be generated via crossover, set in newGeneration as num_individuals minus the number of duplicates
//          rng - random number generator that is used to randomly select parents and is passed into mutate
//          currentAnneal - this is the current anneal status, passed into mutate
//          generation - this is the current generation, passed into mutate
//          cConstants - the cuda constants
// Outputs: This function will fill the newChildren array up to numNewChildren with generated chilren, ready to be simulated
void generateChildrenFromCrossover(std::vector<Adult> & parents, Child* newChildren, const int & numNewChildren, std::mt19937_64 & rng, const double & currentAnneal, const int & generation, const cudaConstants* cConstants);

// This function will use duplicate adults to generate children using heavy mutation
//      It will help ensure that genetic diversity is conserved by not creating children from the crossover between two duplicate parents
// Inputs:  duplicates - a vector of adults that that will used to generate mutated children
//          newChildren - the array that the generated children will be appended to
//          startingIndex - lets the function know where within newChildren to start inserting mutated children
//          rng - the random number generator, will be passed into the mutate function
//          currentAnneal - the mutation anneal, will be passed into the mutate function 
//          generation - the current generation, will be passed into the mutate function
//          cConstants - the cuda constants
// Outputs: Duplicate adults' starting paramters will have been copied into children in the newChildren array and mutated, filling up the rest of newChildren
// NOTE: this function assumes (but isn't dependent on) that generateNewChildrenFromCrossover has been called previously
void generateChildrenFromMutation(std::vector<Adult> & duplicates, Child* newChildren, const int & startingIndex, std::mt19937_64 & rng, const double& currentAnneal, const int & generation, const cudaConstants* cConstants);

// Converts children previously simulated by callRK into adults
//        It calculates the pos and speed diffs of the children and inserts the new children into the submitted adult vector
// Input: newAdults - the vector of adults the converted adults will be inserted into
//        newChildren - the array of simulated children that will be tranformed into adults
//        cConstants - needed for the number of individuals
// Output: newAdults will be filled with the children that are converted into adults 
//this happens after teh first generation is created and then after every new child generation is made 
void convertToAdults(std::vector<Adult> & newAdults, Child* newChildren, const cudaConstants* cConstants);

// Creates the first generation of adults by taking in an array of children with randomly generated parameters
// Input: initialChildren - children with randomly generated parameters whose runge kutta values have not yet been computed 
//        oldAdults - a vector that will be filled with the children once they've gone through callRK and are converted into adults
//        cConstants - passed on into callRK and convertToAdults
// Output: oldAdults is full of the children that have been converted into adults - this allows us to have a bunch of random adults to create a new generation
//
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants);

#include "ga_crossover.cpp"

#endif
