#include "..\Config_Constants\constants.h" //so we can use cuda constants
#include "..\Genetic_Algorithm\adult.h" //allows us to use adults
#include <vector>   // allows us to use vectors instead of just arrays

#ifndef ANNEAL_H
#define ANNEAL_H

//----------------------------------------------------------------------------------------------------------------------------
//Input: current anneal and dRate
//Output: an update of the anneal and dRate based on the tolerance and changeInBest
//Function that adjusts the anneal based on the current circumstances
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday, double & dRate, const int & generation, double & previousProgress);

#include "anneal.cpp"

#endif