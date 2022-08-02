#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

#ifndef ANNEAL_H
#define ANNEAL_H

///////////////////////////////////////////////////////////////
// Anneal Functions                                          //
///////////////////////////////////////////////////////////////

// --- This file includes functions used to calculate anneal ---
// TODO: I would put more info here regarding ^^^^^^^^^^^^^^

//----------------------------------------------------------------------------------------------------------------------------
//This function will calculate a new anneal based on the current progress of the 0th (best) individual
//  NOTE: this function assumes that the inputted vector is sorted in rank-distance order
//Input:  oldAdults - the array of adults from which the best adult will be referenced from
//        cCOnstants - the cuda constants
//        currentAnneal - the generation's anneal value, which will be modified by the function
//        generation - the current generation
// Output: an update of the anneal based on the best current progress
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, const int & generation);

#include "anneal.cpp"

#endif