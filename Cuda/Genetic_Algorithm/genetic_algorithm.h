#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <iostream>

#include "..\Config_Constants\constants.h"
#include "..\Genetic_Algorithm\adult.h"
#include "..\Genetic_Algorithm\child.h"
#include "..\Genetic_Algorithm\ga_crossover.h"
#include "..\Genetic_Algorithm\sort.h"

// newGeneration will generate a mask and use random adults from oldAdults to generate new children
//      It will then simulate the newChildren.
//      Finally, it will convert the new children into adults and insert them into the newAdults vector
//      Note: the number of new children generated should ideally equal the number of threads on the gpu
//          this is set in the config file
// Inputs:  oldAdults - a vector of adults that will be used to pull random parents to generate children (should be at least survivor_size, preferably num_individuals size)
//          newAdults - a vector of adults that will be filled by children that have been simulated and converted into adults (will be num_individuals size)
//          annealing - a double that will be passed on to other functions; it determines how much a child's parameters are mutated
//          generation - the generation that the algorithm is on; this will be passed on to other functions and is used for reporting
//          rng - a random number generator that is both passed on to other functions and used to pull random parents
//          cConstants - the config constants; it is passed on to other functions and is used for constructing children and mutating their parameters
// Outputs: newAdults will be filled with newly generated adults; it will be ready to be sorted and compared to other generations
// NOTE: This function is called at the beginning of each generation within optimize() 
void newGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants);

// Copies unique values from oldAdults into parents and any adults with the same posDiff and speedDiff as a parent gets copied into duplicates
// Input: oldAdults - the potential parents from which we are selecting the best survivor_count parents and duplicates
//        parents - the unique individuals or the first instance of a duplicate
//        duplicates - individuals with the same posDiff or speedDiff as a parent
//        generation - if it is the first generation, nothing will be set with duplicate error status yet, so everything is a parent
//        cConstants - so we can access survivor_count
// Output: parents and duplicates are filled with values copied over from oldAdults
void separateDuplicates(std::vector<Adult> & oldAdults, std::vector<Adult> & parents, std::vector<Adult> & duplicates, const int & generation, const cudaConstants* cConstants);

// Takes parents and duplicates and uses them to fill newChildren with the num_individuals children
// Input: parents - the unique individuals or the first instance of a duplicate
//        duplicates - individuals with the same posDiff or speedDiff as a parent
//        newChildren - the array of simulated children that will be transformed into adults
//        annealing - a double that will be passed on to other functions; it determines how much a child's parameters are mutated
//        generation - the generation that the algorithm is on; this will be passed on to other functions and is used for reporting
//        rng - a random number generator that is both passed on to other functions and used to pull random parents
//        cConstants - the config constants; it is passed on to other functions and is used for constructing children and mutating their parameters
// Output: newChildren is filled with 
void makeChildren(std::vector<Adult> & parents, std::vector<Adult> & duplicates, Child * newChildren, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants);

// Converts children previously simulated by callRK into adults
//        It calculates the pos and speed diffs of the children and inserts the new children into the submitted adult vector
// Input: newAdults - the vector of adults the converted adults will be inserted into
//        newChildren - the array of simulated children that will be transformed into adults
//        cConstants - needed for the number of individuals
// Output: newAdults will be filled with the children that are converted into adults 
//this happens after teh first generation is created and then after every new child generation is made 
void convertToAdults(std::vector<Adult> & newAdults, Child* newChildren, const cudaConstants* cConstants);

// Creates the first generation of individuals within inputParameters
// Input: oldAdults - an empty vector that will be filled with the first set of parents (randomly generated or from a file)
//        cConstants - allows us to access num_individuals and other important information
//        rng - a random number generator so we can get random parameters for the children that will become the first generation of adults
//        generation - child constructor needs a creation generation for its birthday
// Output: oldAdults is filled with the first generation of adults with random parameters
//         calls firstGeneration to turn the randomly generated children into Adults  
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng, const int & generation);

// Creates the first generation of adults by taking in an array of children with randomly generated parameters
// Input: initialChildren - children with randomly generated parameters whose runge kutta values have not yet been computed (from createFirstGeneration)
//        oldAdults - a vector that will be filled with the children once they've gone through callRK and are converted into adults
//        cConstants - passed on into callRK and convertToAdults
// Output: oldAdults is full of the random children that have been converted into adults -> this will be num_individuals adults
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants);

// Sorts the adults from oldAdults and newAdults into allAdults
// Input: allAdults - holds the oldAdults and newAdults from the last generation initially (this is where oldAdults was selected from). 
//                    After this function, it is filled with all the potential parents from this generation
//                    and their offspring, and it is in rankDistanceSort order 
//                    All the bad individuals have been removed from the population, so this is not necessarily 2*num_individuals size
//        oldAdults - passed in with the current parents who created the most recent generation of children
//                    After, this holds the best num_individuals adults from allAdults in rankDistanceSort order (the new parents)
//        newAdults - initially is full of the most recently generated set of children (num_individuals size)
//                    After the function is complete it is empty and ready to hold more children
//        numNans - this function tallies the nans (and other errors) that occur in this population of adults
//        cConstants - allows us to access num_individuals
//        generation - this is needed for eliminateBadAdults
// Output: oldAdults is filled with the adults that are potential parents for the next generation
//TODO: more info about what is happening to allAdults / newAdults / oldAdults
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants, const int & generation);

//function that sorts the adults in newAdults and oldAdults and puts them in allAdults while excluding unwanted adults
//the unwanted adults are those with nans, uneven genes (really good speedDiff, but worst posDiff), and adults that are older than 50 generations
//Input: allAdults - holds the adults we want to keep
//       newAdults - the adults that were just created and need to be checked
//       numNans - this function tallies the nans (and other errors) that occur in this population of adults
//       cConstants - needed for the pos_threshold and speed_threshold
//       generation - to check the current generation to eliminate the old adults
//Output: allAdults is filled with the adults we want 
//this is called in preparePotentialParents()
void eliminateBadAdults(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants, const int & generation);

#include "genetic_algorithm.cpp"

#endif