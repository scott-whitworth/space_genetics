#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

#ifndef ANNEAL_H
#define ANNEAL_H

//----------------------------------------------------------------------------------------------------------------------------
//Input: current anneal and dRate
//Output: an update of the anneal and dRate based on the tolerance and changeInBest
//Function that adjusts the anneal based on the current circumstances
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, int & generation, const double & posTolerance, double & dRate);

//----------------------------------------------------------------------------------------------------------------------------
// Used to identify if the best adult has changed
// Input: previousBestPosDiff - the posDiff of the previous best adult
//        currentBest - the best Adult currently identified
//        distingushRate - utility variable used to help determine if the current best and previous best is different, will allow for comparison, regardless of how small the posDiffs are
// Output: Returns true if there was a change in the best individual
bool changeInBest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate);

//----------------------------------------------------------------------------------------------------------------------------
//This function will calculate the cost based on the mission type and the current best individual
//      This function should only be used within changeAnneal 
//Inputs: oldAdults - A vector of adults that will pull the individual whose final position will be used to calculate the cost
//              NOTE: This function assumes that oldAdults is already sorted by rankDistance and will pull the first adult from the vector
//Output: A cost double that will be used to change anneal              
double calculateCost(const std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

#include "anneal.cpp"

#endif