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

// --- This file includes functions that deal with creating crossover masks and generating children based on the masks ---

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

//Sets the whole mask to be set to AVG_RATIO
// ** currently not in use ** (Tests should be done pertaining to its use)
// Input: mask - pointer integer array of length OPTIM_VARS
// Output: mask is set to contain all AVG_RATIO values
void crossOver_averageRatio(std::vector<int> & mask);

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
//        mutationScale - scalar for the mutation intensity, allows for control on the intensity of mutations between children and duplicates
//        mutation_chance - the chance for each gene to mutate, different if its a duplicate or not
// Output: Returns rkParameter object that is the mutated version of p1
// Called by generateNewIndividual
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, const double & annealing, const cudaConstants* cConstants, const double & mutationScale, const double & mutation_chance);

// Utility function for mutate() to get a random double with high resolution
// Input: max - the absolute value of the min and max(min = -max) of the range
// Output: A double value that is between -max and +max
double getRand(double max, std::mt19937_64 & rng);

///////////////////////////////////////////////////////////////////////////////

// Creates a new rkParameters individual by combining properties of two parent Individuals using a crossover mask
// Input: two rkParameter from individuals (p1 and p2) - source of genes for new individual
//        mask - Contains maskValue values and length of OPTIM_VARS,
//               determines how the two parent properties are merged into creating the new individual
//        cConstants, annealing, rng - passed through to mutate()
// Output: Returns rkParameter object that is new individual
// Called from generateChildrenPair, calls mutate
rkParameters<double> generateNewChild(const rkParameters<double> & p1, const rkParameters<double> & p2, const std::vector<int> & mask, const cudaConstants * cConstants, const double & annealing, std::mt19937_64 & rng);

// Method that creates a pair of new Children from a pair of parent Adults and a mask
// Input:  parent1 - the first parent that the children will draw parameters from
//         parent2 - the second parent that the children will draw parameters from
//         newChildren - a pointer to the newChildren array that is used to add the newly generated children to
//         mask - pointer vector of maskValues used to decide on which property from which parent is acquired (or average of the two)
//         annealing - double variable passed onto mutateNewIndividual
//         rng - random number generator passed on to generateNewChild
//         numNewChildren - Used within newGeneration to keep track of how many children have been generated; added to whenever a new child is generated
//         generation - passed on to GenerateNewChild... Will eventually be used in reporting mutations (currently not being passed to mutate, just being used in child for birthday)
//         cConstants - passed on to generateNewChildren... Used for calculating mutations and constructing children
// Output: The newChildren array will contain two newly generated children
//         mask is flipped in polarity (refer to flipMask method) 
//         numNewChildren is incremented by +2
// Called by newGeneration() each time a new mask is generated 
void generateChildrenPair(const Adult & parent1, const Adult & parent2, Child * newChildren, const int & childrenToGenerate, std::vector<int> & mask, std::mt19937_64 & rng, int & numNewChildren, const int & generation, const cudaConstants* cConstants);

///////////////////////////////////////////////////////////////////////////////

// This function will generate a certain number of children via the ga_crossover method
//      It will ensure that duplicate adults within the survivor pool do not create children
// Inputs:  parents - this adult vector holds the non-duplicate parent list
//          newChildren - this is the array that generated children will be inserted into
//          childrenToGenerate - tracks how many children needs to be generated via crossover, set in newGeneration as num_individuals minus the number of duplicates
//          rng - random number generator that is used to randomly select parents and is passed into mutate
//          currentAnneal - this is the current anneal status, passed into mutate
//          generation - this is the current generation, passed into child
//          cConstants - the cuda constants
// Outputs: This function will fill the newChildren array up to childrenToGenerate with generated children, ready to be simulated
void generateChildrenFromCrossover(std::vector<Adult> & parents, Child* newChildren, const int & childrenToGenerate, std::mt19937_64 & rng, const int & generation, const cudaConstants* cConstants);


#include "ga_crossover.cpp"

#endif
