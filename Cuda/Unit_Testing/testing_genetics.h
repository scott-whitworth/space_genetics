#define UNITTEST
#include "../Planet_calculations/planetInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/adult.h" //to access Adults
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // ONLY in here so I don't get errors when compiling since I am including ga_crossover.h which needs this
#include "../Genetic_Algorithm/ga_crossover.h" //to access the genetics crossover algorithms being unit tested here
#include "../Genetic_Algorithm/genetic_algorithm.h" //for the rest of genetics algorithms being tested
#include "../Genetic_Algorithm/sort.h" //allows me to access giveRank and giveDistance directly
#include <math.h>  //allows us to access isnan
#include <random> //for rng
#include <time.h>

//The number of children currently generated by 1 pair of adults by crossover
const int crossoverChildrenCount = 6;
//the offset of the 50/50 average in an array of children generated by crossover 
//important because this seems like the best way for determining an individuals' parents
const int fullAvgOffset = 2;

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

// Takes in the adults from firstParentsTest and checks them against the expected results
// Input: theResults - the Adults generated in firstParentsTest
// Output: returns true if the results match our expectations or false if they do not
bool checkParentsTest(std::vector<Adult>& theResults);

// A function that determines if masks are being generated correctly
// Input: rng - a psuedorandom number generator used by to create masks
//        printMask - whether or not you want the numeric values contained in the mask printed
// Output: Outputs true if the masks were all created as expected or false if they weren't
//         If printMask is true, it will also print all the numeric values contained in each mask
bool createMasks(std::mt19937_64& rng, bool printMask);

//makes children using the different method -> not to populate a generation or anything, just to check each child creation method works as expected
//starts with two known parents with know 
// Input: rng - a random number for the mask functions
//        utcConstants - allows you to access mutation rate, thruster type, and is used by create children pair
//        printThings - this will print some more of the steps so you can independently verify the results in case the unit test missed something
// Output: Makes 6 different children (with thruster on and off) and with no mutations 
//         It will return true if these children were generated correctly or false if they were not
bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, const cudaConstants * utcConstants, bool printThings);

// A function that says if the tripTimes, alphas, betas, and zetas are reasonable for the different children
// Input: c1 & c2 - the children generated by a certain method in makeChildrenWithDifferentMethods
//        mask - the mask that was used to make the children
//        parentValues - a vector full of the parents' values for alpha, beta, zeta, and tripTime
//        whichMethod - 1 stands for crossOver_wholeRandom, 2 stands for crossOver_average, 3 stands for crossOver_bundleVars
//        utcConstants - so the function can access thruster type
//        printThings - if you want to see that the code ran successfully for each test, set printThings to true
// Output: True if the output of the code is as expected, false if it is not
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double* parentsValues, int whichMethod, cudaConstants * utcConstants, bool printThings);

// Helper for checkReasonability
// Gets a child's value for alpha, beta, zeta, or triptime, by taking in the offset and a child
// Input: correspondingOffset - ALPHA_OFFSET, BETA_OFFSET, ZETA_OFFEST, or TRIPTIME_OFFSET 
//                              allows you to access a child's alpha, beta, zeta, or triptime
//        aChild - the Child you want information from
// Output: the value of the parameter corresponding to that offset in the child
double getParamStuff(const int & correspondingOffset, const Child& aChild);

// Fills allAdults with 20 adults with determined posDiffs, speedDiffs, and tripTimes 
// Input: printThings - TRUE-> prints the adults after they have been rankDistance sorted
//                     FALSE-> skips printing
//       allAdults - vector of adults overwritten by the adults contained in this function
//       utcConstants - used to access num_individuals (which this function increases to 20)
// NOTE: None of these Adults went through callRK and some of their tripTimes fall outside the acceptable range
//       These values were picked largely as they make the calculations easier and the final version using callRK
//       is tested in firstFullGen
// Output: num_individuals is 20 and allAdults is filled with the 20 individuals that were created in this function
void twentyAdultsPosAndSpeedDiffMade(bool printThings, std::vector<Adult>& allAdults, cudaConstants* utcConstants);

// Checks if clones are being separated from original/unique values properly
// Input: printThings - TRUE-> displays information on the values being held in duplicates
//                      FALSE-> no cout statements are printed unless there is an issue with the code
//        utcConstansts - used to access and change survivor_count, num_individuals, and to pass into the genetic algorithms
// Output: The parents and duplictaes vectors are compared to the expected versions of these vectors 
//         (these vectors were determined by looking at the order things ended up in after rankDistanceSort)
//         Returns TRUE if the values match their expected order
//         FALSE if not
bool verifyProperCloneSeparation(bool printThings, cudaConstants* utcConstants);

// Verifies that generateChildrenFromCrossover works as expected by using different parents
//      and different proportions of these parents to generate children 
//      these children have their unmutated tripTimes compared to their potential parents' to see if they match expectations
// Input: rng - a psuedorandom number generator used by to create masks
//        printThings - TRUE -> prints all the Adults
//                      FALSE -> only prints error messages and whether it passed or failed
//        utcConstants - used to ensure num_individuals and anneal_initial are the desired values
//                       also used to update survivor_count and to pass into genetic algorithms
// Output: TRUE - if all the children generated seem to match the values we expected them to hit
//         FALSE - if any of the children generated are unusual or fail to meet our expectations for them
bool verifyChildrenFromCrossover(std::mt19937_64& rng, bool printThings, cudaConstants* utcConstants);

// Compares the results from verifyChildrenFromCrossover with the expected values for those Children
// Tests a variety of cases with different survivor population sizes and different proportions of unique individuals and duplicates
// Input: endSpot - if the number of children generated using this method is not num_individuals
//                  this Child can be used to verify that spots at the end of the array are not being overwritten
//                  It is the Child that is past the last spot that should be filled by generateChildrenFromCrossover()
//        numChildren - the number of children that were supposed to be generated using generateChildrenFromCrossover()
//        childrenGenerated - the array holding all the children generated 
//        utcConstants - allows us to access num_individuals
//        parents - the vector holding the parents that generated numChildren the children in childrenGenerated array.
//        printThings - TRUE -> prints the parents' tripTimes and the tripTimes of their children
//                      FALSE -> nothing is printed from this function
// Output: TRUE - if all the children match the values we expect them to and take up the proper amount of space in the vector
//                The way we verify our expectations for the children is by verifying that they had two unique parents
//                this is generally done by comparing the Child that should come from the crossOver_average mask's tripTime 
//                to the average tripTime of two potential parents -> the averages of the tripTime values should be unique so 
//                it is possible to verify that a Child was generated from the two parents we think it came from.
//                If not enough children have been generated to look at one generated from the crossOver_average mask, 
//                the code will attempt to determine the parents by looking at the children generated from crossOver_wholeRandom 
//                (this method may be inconclusive, so while it will return true, it also prints a message that it is inconclusive)
//         FALSE - if any of the children do not match their expected values or if any of the elements of the array that are
//                 supposed to remain constant change (if children are overwritten when they shouldn't be)
bool cfcAnswersMatchExpectations(const int & numChildren, const Child* childrenGenerated, const cudaConstants* utcConstants, std::vector<Adult>& parents, bool printThings, int duplicateCount);

// Ensures generateChildrenFromMutation is working as expected
// NOTE: while this function does pass make a set of children using this method when the mutation chance is 0 (both for the default and duplicates),
//           it does NOT actually ensure that none of the parameters were mutated, just that not all of them were
// Input: rng - used to create the mask
//        printThings - determines if extra cout statements (including the status of all the parameters in duplicates and their children)
//                      should be printed or not
//        utcConstants - allows us to access and change anneal, mutation chances, survivor_count, etc
//                       also needed for a lot of the genetic algorithms
// Output: TRUE - all the tests passed successfully and returned the results we expected
//         FALSE - at least one of the tests returned an unexpected value
bool verifyChildrenFromMutation(std::mt19937_64& rng, bool printThings, cudaConstants* utcConstants);

// Used to check the results from verifyChildrenFromMutation
// Input: duplicates - all the duplicates that have produced children
//        children - an array of children, where the last children were generated from the duplicates using generateChildrenFromMutation
//        utcConstants - allows us to access num_individuals and those sorts of things
//        printThings - TRUE -> the elements of startParams of the duplicate and its child that may be mutated are printed
//                      FALSE -> no visual output from this function, you rely more heavily on it being coded correctly
// Output: TRUE - every child has been mutated at least somewhat from the duplicate is was created from
//         FALSE - at least one Child is a duplicate of the duplicate it was created from
bool properMutationResults(std::vector<Adult>& duplicates, Child* children, const cudaConstants* utcConstants, bool printThings);

// Made to ensure they work well enough to use to use in properMutationResults
// Turns out that the compare functions will not work very well (even after I editted them so they could handle nans)
//      it does not do well with small changes in decimals and pretty much all the mutated numbers are decimals that may have small changes
//      so they will not be used in properMutationResults  
// Because these compare functions do not appear to be used anywhere in our code and will not be used here, they are not fully tested here 
// Input: printThings: TRUE - prints a the actual values for compare in a couple places
//                     FALSE - only prints whether or not it worked as expected in a given scenario
// Output: FALSE - these will not work well in properMutationResults
//         TRUE - these would work well in properMutationResults 
bool testingCompareFunctions(bool printThings);

// Uses UTCopyOfNewGen to generate a vector of newAdults 
// These Adults are then converted back into children so that the functions used to verify that mutation and crossover were working previously can be used again
// Input: printThings - TRUE -> prints parents and the children they produced (a lot of stuff)
//                      FALSE -> does not print anything except maybe an error message or so if applicable
//        rng - needed for UTCopyOfNewGen to make masks and do mutations, etc
//        utcConstants - allows us to access and edit num_individuals, survivor_count, mutation rates, etc.
// Output: TRUE -> the newAdults (as Children) were generated in such a way that the verification methods say they were correctly generated
//         FALSE-> newAdults (as Children) had an unexpected error
bool testingCopyOfNewGen(bool printThings, std::mt19937_64 & rng, cudaConstants* utcConstants);

// This function was creared to convert Adults back into children so they could be run through the old verification methods
// Input: newAdult - the Adults that need to be converted to Children
//        newChildren - initially an empty array, but wull be filled with all the adults that have been converted into children
//        utcConstants - to access survivor_count and stuff like that
// Output: newChildren is full of the newAdults, but as Children
void convertBackToChildren(std::vector<Adult>& newAdult, Child* newChildren, const cudaConstants* utcConstants);

//copied directly from ga_crossover.cpp, just doesn't mess with the whole callRK thing
// Makes a new generation of Adults from a vector of old Adults
// Inputs:  oldAdults - a vector of adults that will be used to pull random parents to generate children (should be at least survivor_size, preferably num_individuals size)
//          newAdults - a vector of adults that will be filled by children that have been simulated and converted into adults (will be num_individuals size)
//          annealing - a double that will be passed on to other functions; it determines how much a child's parameters are mutated
//          generation - the generation that the algorithm is on; this will be passed on to other functions and is used for reporting
//          rng - a random number generator that is both passed on to other functions and used to pull random parents
//          utcConstants - the config constants; it is passed on to other functions and is used for constructing children and mutating their parameters
// Outputs: newAdults will be filled with newly generated adults; it will be ready to be sorted and compared to other generations
int UTCopyOfNewGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* utcConstants);

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
bool firstFullGen(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings, GPUMem & gpuValues);

// A function that is used to verify that firstFullGen is working correctly
// Attempts to determine the parents of a set of children by looking at the child made from a whole average
// Input: youngGen - the new Adults that have just been generated
//        possParents - the vector of Adults that could potentially be the parents of the Adults
//        utcConstants - the CUDA constants so we can access num_individuals, num_survivors, and those sorts of things
//        printThings - denotes whether the user wants intermediate steps printed or only a pass/fail message
// Output: TRUE - if the code is able to identify two unique parents for every six Adults in youngGen, using the tripTime of the 3rd Adult
//         FALSE - if the code cannot identify two unique parents for a set of Adults
bool verifyFullGen(std::vector<Adult>& youngGen, std::vector<Adult>& possParents, const cudaConstants * utcConstants, bool printThings);

// Basically the same as firstFullGen, except it uses 7 and 5 parents to generate 65 children
// Input: rng - used by the masks to set up which parent a gene is being taken from
//        youngGen - the last generation created by possParents that should be overwritten by the new generation
//        possParents - the set of potential parents for a generation
//        utcConstants - used for a variety of things from num_individuals to being used by callRK and stolenGiveRank
//        printThings - if print things is true, it will print the parents' tripTimes, highlight the parents that were 
//                      selected to make the next generation, and print the tripTimes for the new generation 
// Output: TRUE - using both 7 and 5 parents, sets of 65 children could be generated and could have their parents identified
//         FALSE - it either fails to identify the parents of 65 children when it either has 7 parents (returns false immediately and doesn't check if it will work for 5)
//                 or returns false if it fails when there are 5 survivors
bool makeManyChildren(std::mt19937_64& rng, std::vector<Adult>& youngGen, std::vector<Adult>& possParents, cudaConstants * utcConstants, bool printThings, GPUMem & gpuValues);

//NOT CURRENTLY USING THE BELOW TESTS 

// Sets up the mask and mutation_rate, prints how many were set to true, and prints the values in the mask
// Input: NA
// Output: The number of genes mutated with a mutation rate of 0.75 and the mask indicating which genes will be mutated
bool checkUTMutateMask();

// Unit test version of the mutateMask function from ga_crossover 
//      as of summer 2022 basically the same as the actual function, only it has a cout statement
//      at the end that prints the number of genes that were mutated
// Input: rng - random number generating object, defined in optimize()
//        mutateMask - pointer to a boolean array, length of OPTIM_VARS
//                   - sets genes being mutated to true and others to false
//        mutation_rate - a double value less than 1 that is the chance a gene will be mutated
//                      - called iteratively to mutate more genes
// Output: mutateMask contains false for genes that are not mutating, true for genes that are to be mutated
//      Prints how many genes were mutated
void UTmutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate);

#include "testing_genetics.cpp"