#ifndef SORT_H
#define SORT_H

#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

///////////////////////////////////////////////////////////////
// Sort Functions                                            //
///////////////////////////////////////////////////////////////

// --- This file includes functions that handle assigning rank and distance to adults so they can be sorted ---

//----------------------------------------------------------------------------------------------------------------------------
//Used to give rankings for sorting based on non-dominated sorting method. Used for rendezvous mission
//Assigns suitability rank to all adults.
//Inputs: allAdults - this generation of adults, will be compared to eachother to find rank
//        cConstants - the cuda constants
//Output: the adults in allAdults will have been assigned a rank
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Gives a distance value to each adult. A higher distance indicates that it is more diverse from other adults
// The distance is how different an adult is from the two adults closest to it for each objective function.
// Inputs: allAdults - the vector of adults which will have their distances calculated
//         cConstants - the cuda constants
// Output: Calculates the distance for all individuals in the pool
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Assist function for objectives that sorts an adult vector based on the goal 
//  This will specifically look at the objective's goal and choose a non-rank/rankDistance sort based on it 
//  As of writing this function, it is only used for outputs and calculating distance
// Inputs: adults - the vector of adults that will be sorted based on the objective
//         sortObjective - the objective the sort will be derived from
//         sortSize - how many adults of the passed in vector will be sorted 
//              The sorting range will span from the beginning of adult to adult.begin() + sortSize
// Output: the adult vector will be sorted with the parameter and order based on the goal of the objective
void parameterSort(std::vector<Adult> & adults, const objective& sortObjective, const int & sortSize);

#include "sort.cpp"

#endif