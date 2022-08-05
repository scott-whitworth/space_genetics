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

///////////////////////////////////////////////////////////////
// Genetic Algorithim Functions                              //
///////////////////////////////////////////////////////////////

// For more information on this section of code, please refer to the geneticAlgorithm.cpp/h section of Code_Overview_CS_readme.md

// --- This file includes functions which handle generation management, such as selecting parents and generation creation ---

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
//          gpuValues - allows callRK to access the the GPU memory that has been allocated and the set values for Mars and the cConstants
// Outputs: newAdults will be filled with newly generated adults; it will be ready to be sorted and compared to other generations
// NOTE: This function is called at the beginning of each generation within optimize() 
void newGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants, GPUMem & gpuValues);

// Copies unique values from oldAdults into parents (size N/4) and any adults with the same a non-duplicate error status are included
// Input: oldAdults - the potential parents from which we are selecting the best survivor_count parents
//        parents - the best individuals
//        generation - if it is the first generation, everything is a parent
//        cConstants - so we can access survivor_count
// Output: parents are filled with values copied over from oldAdults
void fillParents(std::vector<Adult> & oldAdults, std::vector<Adult> & parents, const int & generation, const cudaConstants* cConstants);

// Takes parents and uses them to fill newChildren with the num_individuals children
// Input: parents - the best N/4 individuals
//        newChildren - the array of simulated children that will be transformed into adults
//        annealing - a double that will be passed on to other functions; it determines how much a child's parameters are mutated
//        generation - the generation that the algorithm is on; this will be passed on to other functions and is used for reporting
//        rng - a random number generator that is both passed on to other functions and used to pull random parents
//        cConstants - the config constants; it is passed on to other functions and is used for constructing children and mutating their parameters
// Output: newChildren is filled with 
void makeChildren(std::vector<Adult> & parents, Child * newChildren, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants);

// Converts children previously simulated by callRK into adults
//        It calculates the pos and speed diffs of the children and inserts the new children into the submitted adult vector
// Input: newAdults - the vector of adults the converted adults will be inserted into
//                    NOTE: newAdults is cleared in the beginning of the function, it is not needed for this to be empty before
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
//        gpuValues - needs to be passed into firstGeneration so callRK can use it to access GPU memory that has already been allocated
// Output: oldAdults is filled with the first generation of adults with random parameters
//         calls firstGeneration to turn the randomly generated children into Adults  
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng, const int & generation, GPUMem & gpuValues);

// Creates the first generation of adults by taking in an array of children with randomly generated parameters
// Input: initialChildren - children with randomly generated parameters whose runge kutta values have not yet been computed (from createFirstGeneration)
//        oldAdults - a vector that will be filled with the children once they've gone through callRK and are converted into adults
//        cConstants - passed on into callRK and convertToAdults
//        gpuValues - allows callRK to access the the GPU memory that has been allocated and the set values for Mars and the cConstants
// Output: oldAdults is full of the random children that have been converted into adults -> this will be num_individuals adults
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants, GPUMem & gpuValues);

// Sorts the adults from oldAdults and newAdults into allAdults
// Input: allAdults - holds the oldAdults and newAdults from the last generation initially (this is where oldAdults was selected from). 
//                    After this function, it is filled with all the potential parents from this generation
//                    and their offspring, and it is in rankDistanceSort order 
//                    All the bad individuals have been removed from the population, so this is not necessarily 2*num_individuals size
//        oldAdults - passed in with the current parents who created the most recent generation of children
//                    After, this holds the best num_individuals adults from allAdults in rankDistanceSort order (the new parents)
//        newAdults - initially is full of the most recently generated set of children (num_individuals size)
//                    After the function is complete it is empty and ready to hold more children
//        numErrors - this function tallies the nans (and other errors) that occur in this population of adults
//        cConstants - allows us to access num_individuals
//        generation - this is needed for eliminateBadAdults
//        currentAnneal - this is used in findDuplicate to change the tolerance for finding a duplicate
// Output: oldAdults is filled with the top N number of adults, sorted by rankDistance
//              The adults that are potential parents for the next generation
//         allAdults is filled with 2N - number of duplicates + errors adults from the combination of the inputted old/newAdults vectors
//              It is sorted by rankDistance
//         newAdults remains the same (size of N and unsorted)
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numErrors, int& duplicateNum, const cudaConstants* cConstants, const int & generation, const double& currentAnneal, int& marsErrors);

//function that decides whether or not to add an adult to allAdults based on its errorStatus
//Input: adultPool - This is either new or oldAdults and an adult from the vector will be checked for errorStatus and added to allAdults
//       allAdults - holds the adults we want to keep
//       index - index of the Adult being checked
//       numErrors - this function tallies the nans (and other errors) that occur in this population of adults
//       duplicateNum - counts the duplicates in both new and oldAdults
//Output: allAdults is filled with the adults we want and counts are updated
//this is called in preparePotentialParents()
void addToAllAdults(std::vector<Adult> & adultPool, std::vector<Adult> & allAdults, const int & index, int& numErrors, int& duplicateNum, int& marsErrors);

#include "genetic_algorithm.cpp"

#endif