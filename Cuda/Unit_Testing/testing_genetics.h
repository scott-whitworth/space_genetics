#define UNITTEST
#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/adult.h" //to access Adults
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // ONLY in here so I don't get errors when compiling since I am including ga_crossover.h which needs this
#include "../Genetic_Algorithm/ga_crossover.h" //to access the genetics algorithms being unit tested here
#include <random> //for rng
#include <time.h>

//const int genSize = 10; //a now unnecessary variable basically representing num_individuals (used in wrongWayToRank which is no longer in use)

// Makes a set of cudaConstants that so ga_crossover functions could be used directly and calls all the unit test functions
// Input: printThings
// Output: Makes a default set of cudaConstants
bool runGeneticsUnitTests(bool printThings);

// Returns true if the first generation of parents is generated from a set of "random" children
// These parents cannot actually be used to create the next generation, but allow for checking giveRank, giveDistance, and rankDistanceSort...
//       ...when you have a set of Adults with known posDiffs and speedDiffs without needing to use rkParameters
// NOTE: I could not include the actual giveRank and giveDistance from optimization.cu since that also has an int main so I copied them over
//      as stolenGiveRank and stolenGiveDistance (their headers specify the last time they were updated)
// Input: utcConstants - used to access num_individuals and used for giveRank and giveDistance
//        printThings - this will print some more of the steps so you can independently verify the results in case the unit test missed something
// Output: returns true if this generation can be correctly generated and sorted or false if not
bool firstParentsTest(const cudaConstants * utcConstants, bool printThings);

// A function that determines if masks are being generated correctly
// Input: rng - a psuedorandom number generator used by to create masks
//        printMask - whether or not you want the numeric values contained in the mask printed
// Output: Outputs true if the masks were all created as expected or false if they weren't
//         If printMask is true, it will also print all the numeric values contained in each mask
bool createMasks(std::mt19937_64& rng, bool printMask);

// Takes in the adults from firstParentsTest and checks them against the expected results
// Input: theResults - the Adults generated in firstParentsTest
// Output: returns true if the results match our expectations or false if they do not
bool checkParentsTest(std::vector<Adult>& theResults);

//makes children using the different method -> not to populate a generation or anything, just to check each child creation method works as expected
//starts with two known parents with know 
// Input: rng - a random number for the mask functions
//        utcConstants - allows you to access mutation rate, thruster type, and is used by create children pair
//        printThings - this will print some more of the steps so you can independently verify the results in case the unit test missed something
// Output: Makes 6 different children (with thruster on and off) and with no mutations 
//         It will return true if these children were generated correctly or false if they were not
bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings);

// A function that says if the tripTimes, alphas, betas, and zetas are reasonable for the different children
// Input: c1 & c2 - the children generated by a certain method in makeChildrenWithDifferentMethods
//        mask - the mask that was used to make the children
//        parentValues - a vector full of the parents' values for alpha, beta, zeta, and tripTime
//        whichMethod - 1 stands for crossOver_wholeRandom, 2 stands for crossOver_average, 3 stands for crossOver_bundleVars
//        utcConstants - so the function can access thruster type
//        printThings - if you want to see that the code ran successfully for each test, set printThings to true
// Output: True if the output of the code is as expected, false if it is not
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double parentsValues[], int whichMethod, cudaConstants * utcConstants, bool printThings);

// Helper for checkReasonability
// Gets a child's value for alpha, beta, zeta, or triptime, by taking in the offset and a child
// Input: correspondingOffset - ALPHA_OFFSET, BETA_OFFSET, ZETA_OFFEST, or TRIPTIME_OFFSET 
//                              allows you to access a child's alpha, beta, zeta, or triptime
//        aChild - the Child you want information from
// Output: the value of the parameter corresponding to that offset in the child
double getParamStuff(const int correspondingOffset, const Child& aChild);

// Makes a set of parents and creates children from these parents -> uses callRK
// Input: rng - used by the masks to set up which parent a gene is being taken from
//        utcConstants - used for a variety of things from num_individuals to being used by callRK and stolenGiveRank
//        printThings - if print things is true, it will print the parents' tripTimes, highlight the parents that were 
//                      selected to make the next generation, and print the tripTimes for the new generation
//        NOTE: tripTime was selected as it is the easiest thing to see a change in and to compare with other individuals
//        NOTE: the name of the vector "parents" may be misleading - these are really the potential parents and only 
//              survivor_count of them actually produce the next generation
// Output: Returns true if the results seem to be as expected and false if there is an issue with them
//         If printThings is true, it will print some things to the terminal, but even if it is false, it will print any error messages
bool firstFullGen(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings);

// A function that is used to verify that firstFullGen is working correctly
bool verifyFullGen(std::vector<Adult>& youngGen, std::vector<Adult>& possParents, const cudaConstants * utcConstants, bool printThings);

bool makeManyChildren(std::mt19937_64& rng, std::vector<Adult>& youngGen, std::vector<Adult>& possParents, cudaConstants * utcConstants, bool printThings);

//sets up the mask and mutation_rate and prints how many were set to true
bool checkUTMutateMask();

//unit test version of the mutateMask function from ga_crossover
void UTmutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate);

//took this directly from optimization.cu on 6/14/22 near the end of the day - will become out of date if changes made to version in optimization.cu
void stolenGiveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants);
//took this directly from optimization.cu on 6/14/22 near the end of the day - will become out of date if changes made to version in optimization.cu
void stolenGiveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants);


//giveRank was broken so I threw together an inefficient giveRank type function 
//Not sure it does the same thing, but works well enough for me to be able to verify if the results seems reasonable
//Calculates rank based on net dominations- different than actual giveRank -> makes more ranks
//void wrongWayToRank(std::vector<Adult> & newAdults);



#include "../Unit_Testing/testing_genetics.cpp"