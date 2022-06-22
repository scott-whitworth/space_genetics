#ifndef SORT_H
#define SORT_H

#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

//----------------------------------------------------------------------------------------------------------------------------
//Used to give rankings for sorting based on non-dominated sorting method. Used for rendezvous mission
//Assigns suitability rank to all adults.
//Input: pool - this generation of adults, defined/initialized in optimimize
//       cConstants
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Gives a distance value to each adult. A higher distance indicates that it is more diverse from other adults
// The distance is how different an adult is from the two adults closest to it for each objective function.
// Input: pool - the full pool of 2N adults excluding NaNs
//        poolSize - the poolSize of all the adults excluding NaNs
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

#include "sort.cpp"

#endif