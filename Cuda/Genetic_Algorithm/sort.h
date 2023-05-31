#ifndef SORT_H
#define SORT_H

#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays
//#include <algorithm> // for std::sort() in sort.cpp (needed for Tesla only)

///////////////////////////////////////////////////////////////
// Sort Functions                                            //
///////////////////////////////////////////////////////////////

// All sorting / distance / rank calculations based off of https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

// --- This file includes functions that handle assigning rank and distance to adults so they can be sorted ---

//----------------------------------------------------------------------------------------------------------------------------
//Used to give rankings for sorting based on non-dominated sorting method.
//Assigns suitability rank to all adults.
//Inputs: allAdults - this generation of adults, will be compared to eachother to find rank
//                    because they are adults, they are guarenteed to have valid speedDiff/posDiff
//        cConstants - the cuda constants for the given run
//Output: the adults in allAdults will have been assigned a rank
//            rank is a property of Adult
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Gives a distance value to each adult. A higher distance indicates that it is probably more diverse from other adults
// The distance is how different an adult is from the two adults closest to it for each objective function.
// Inputs: allAdults - the vector of adults which will have their distances calculated
//                   - these adults will have already been given rank (giveRank)
//         cConstants - the cuda constants
// Output: Calculates the distance for all individuals in the pool
//             distance is stored as a member of Adult
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Assist function for objectives that sorts an adult vector based on the goal 
//  This will specifically look at the objective's goal and choose a non-rank/rankDistance sort based on it 
//  designed to be used by giveDistance
// Inputs: adults - the vector of adults that will be sorted based on the objective
//         sortObjective - the objective that will be sorted  
//         sortSize - how many adults of the passed in vector will be sorted 
//              The sorting range will span from the beginning of adult to adult.begin() + sortSize
//              this is used to only sort the first part of the vector (giveDistance)
// Output: the adult vector will be sorted with the parameter and order based on the goal of the objective
void parameterSort(std::vector<Adult> & adults, const objective& sortObjective, const int & sortSize);

#include "sort.cpp"

#endif