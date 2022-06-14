#define UNITTEST
//#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/adult.h"
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // ONLY in here so I don't get errors when compiling since I am including ga_crossover.h which needs this
#include "../Genetic_Algorithm/ga_crossover.h"
#include <random>

const int genSize = 10;

// Makes a set of cudaConstants that so ga_crossover functions could be used directly and calls all the unit test functions
// Input: NA
// Output: Makes a default set of cudaConstants
bool runGeneticsUnitTests();

// A function that determines if masks are being generated correctly
// Input: rng - a psuedorandom number generator used by to create masks
//        printMask - whether or not you want the numeric values contained in the mask printed
// Output: Outputs true if the masks were all created as expected or false if they weren't
//         If printMask is true, it will also print all the numeric values contained in each mask
bool createMasks(std::mt19937_64& rng, bool printMask);

//returns true if the first generation of parents is generated from a set of "random" children
//these parents cannot actually be used to create the next generation I realized after the fact because their rkParameters and elements are the defaults 
//I made a child constructor that allowed me to set the speedDiff and posDiff so it doesn't correspond to their rkParameters or elements
//The logic for the above decision was that it would make it easier to tell if the individuals were properly becoming adults by rankDistanceSorting them
bool firstParentsTest(const cudaConstants * utcConstants);

//makes children using the different method -> not to populate a generation or anything, just to check each child creation method works as expected
//starts with two known parents with know 
bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, cudaConstants * utcConstants);

//Are the tripTimes, alphas, betas, and zetas reasonable for the different children
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double parentsValues[], int whichMethod);

//Makes a set of parents and creates children from these parents
bool firstFullGen(const cudaConstants * utcConstants);

//TODO: once giveRankTest and giveDistanceTest are complete, use these instead, rather than using my weird versions of these
//giveRank was broken so I threw together an inefficient giveRank type function 
//Not sure it does the same thing, but works well enough for me to be able to verify if the results seems reasonable
void wrongWayToRank(std::vector<Adult> & newAdults);
//copied giveDistance from optimization.cu and updated it as necessary to work with my other testing functions
void sortaGiveDistance(std::vector<Adult> & pool);

//sets up the mask and mutation_rate and prints how many were set to true
bool checkUTMutateMask();

//unit test version of the mutateMask function from ga_crossover
void UTmutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate);


#include "../Unit_Testing/testing_genetics.cpp"