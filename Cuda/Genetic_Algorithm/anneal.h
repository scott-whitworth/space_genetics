#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

#ifndef ANNEAL_H
#define ANNEAL_H

///////////////////////////////////////////////////////////////
// Anneal Functions                                          //
///////////////////////////////////////////////////////////////

// --- This file includes functions used to calculate anneal ---
//The anneal_initial is initialized in the config, a good value for this is 0.02, anything that starts larger than 0.05 is usually too large.
//The currentAnneal is adjusted every generation in this function.
//A lot of testing of different ways to adjust annealing and the rates of convergences led to this method.
//This method was made for the bennu mission but should work for any mission.
//The main concept behind this method is to avoid points where the code gets stuck for 100s of generations.
//The sine function adds some level of variability that other methods don't have.
//The reason for a step method is its the easiest way to adjust the anneal in small steps, as there is no function that behaves the way we need to properly adjust anneal.
//These steps are also the result of many tests, but are based off nothing but intuition, so they could/should be changed in the future. 
//The last function used for the 1e-3-1.0 range is a function that decays quickly within that range of progress. Again, this function could/should be changed at some point.


//----------------------------------------------------------------------------------------------------------------------------
//This function will calculate a new anneal based on the current progress of the 0th (best) individual
//  NOTE: this function assumes that the inputted vector is sorted in rank-distance order
//Input:  oldAdults - the array of adults from which the best adult will be referenced from
//        cCOnstants - the cuda constants
//        currentAnneal - the generation's anneal value, which will be modified by the function
//        generation - the current generation
// Output: an update of the anneal based on the best current progress
double changeIndividualAnneal (double curProgress, const cudaConstants* cConstants, const int & generation);
#include "anneal.cpp"

#endif