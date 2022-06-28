#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

#ifndef ANNEAL_H
#define ANNEAL_H

//----------------------------------------------------------------------------------------------------------------------------
//Input: current anneal and dRate
//Output: an update of the anneal and dRate based on the tolerance and changeInBest
//Function that adjusts the anneal based on the current circumstances
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday, double & dRate);

//----------------------------------------------------------------------------------------------------------------------------
//This function will calculate the cost based on the mission type and the current best individual
//      This function should only be used within changeAnneal 
//Inputs: oldAdults - A vector of adults that will pull the individual whose final position will be used to calculate the cost
//              NOTE: This function assumes that oldAdults is already sorted by rankDistance and will pull the first adult from the vector
//Output: A cost double that will be used to change anneal              
double calculateCost(const std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

#include "anneal.cpp"

#endif