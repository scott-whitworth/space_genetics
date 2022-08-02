#ifndef ADULT_H
#define ADULT_H

#include "../Genetic_Algorithm/child.h"
#include "../Runge_Kutta/rkParameters.h"
#include "../Planet_calculations/planetInfo.h"
#include <vector>
#include <string> //used for the unitTestingRankDistanceStatusPrint method so we can get rank, distance, and status easily during unit testing

///////////////////////////////////////////////////////////////
// Adult Class & Functions                                   //
///////////////////////////////////////////////////////////////

// --- This file includes the definitions for the Adult class and associated functions ---

//NOTE FOR FUTURE UNDERSTANDING: the term individuals within comments is used as a general term for both adults and children.
//                               Individual most often refers to adults.
//                               A good rule of thumb: Adults should never be simulated, Children should never be sorted

// Adult is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
// Adult is inheriting startParams, errorStatus, posDiff, speedDiff, and finalPos from Child
// adult inherits the genes from child to avoid calling RK when we dont need to
struct Adult: public Child {
    //rank measures how good this adult is compared to other adults in it's generation. The lower the rank, the better.
    //     NOTE: ranks are not exclusive to one individual (e.g. multiple individuals can have rank 1); it may be helpful to think of ranks as tiers instead
    //It is set within the giveRank() function in optimization.cu
    //     NOTE: DO NOT try to call rank before using giveRank(), even for oldAdults; oldAdults wouldn't have been compared to the new individuals yet, so if they already have a rank, it is innacurate
    //Rank is used as a basis to sort individuals, specifically within rankSort and rankDistanceSort
    //It is also helpful in determining the current best individual 
    int rank; 

    //distance is how different this adult is from the adults most similar to it.
    //It is calculated using the giveDistance() function in optimize.cu
    //It is used within rankDistanceSort to differentiate between adults of the same rank
    //A large distance can be thought of as a high degree of genetic diversity for the adult
    double distance; 

    //this is changed/checked in newGeneration()
    //meant to make sure there is not more than one instance of the same adult in the next generation


//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    // Input: none
    // Output: an adult with a maximum rank and a distance of -1
    // NOTE: This should not be used, it is here as a guardrail for if we are unintentionally generating adults
    Adult();

    // Constructor to create an adult from a child
    // Input: A child - this has cConstants, a set of rkParameters, and elements in it already, as well as the speedDiff and posDiff
    // Output: a child is turned into an adult
    // NOTE: This will be called by convertToAdults within newGeneration & ga_crossover.cu. Children that have just been simulated will be converted from a child into an adult
    Adult(const Child& c): Child(c), rank(INT_MAX), distance(-1){}

    #ifdef UNITTEST //this should be defined in unit testing

    // Constructor used ONLY FOR UNIT TESTING!!!
    // This constructor uses the default rkParameters, elements, posDiff and speedDiff, the only things that are defined are the rank and distance because that is all rankDistance sort cares about
    // Input: A rank (r), a distance (d), and a STATUS (s) -> status defaults to VALID if undefined
    // Output: Constructs an Adult using the unit testing Child constructor that only contains default and sets its rank and distance so it can be sorted using rankDistanceSort
    // ONLY USE FOR UNIT TESTING!!! SHOULD NEVER BE CALLED ELSEWHERE 
    Adult(int r, double d, ERROR_STATUS s = VALID): Child(s), rank(r), distance(d){}

    // A function used exclusively for unit testing that prints that rank and distance of an adult
    // Output: takes an adult's rank, distance, and status and creates a string holding this information
    // ONLY USE FOR UNIT TESTING!!! Note: This was purposefully given a long name to discourage people from using it
    const std::string unitTestingRankDistanceStatusPrint();

    //getters used exclusively in unit testing
    int getRank();
    double getDistance();
    #endif //unit testing endif

};

//TODO: We need to Refactor the below if we are going to change the optimized parameters (velocity direction)

// Compare two adults by their positional difference values
// input: two adults
// output: returns true if personB has a higher positional difference than personA
bool LowerPosDiff(Adult& personA, Adult& personB);

// Compare two adults by their velocity difference values
// input: two adults
// output: returns true if personB has a higher velocity difference than personA
bool HigherSpeedDiff(const Adult& personA, const Adult& personB);

// Compare two adults by their velocity difference values
// input: two adults
// output: returns true if personB has a lower velocity difference than personA
bool LowerSpeedDiff(const Adult& personA, const Adult& personB);

// Compare two individuals by their spent fuel values, used in standard sort to have array contain lowest fuel spent individual at start
// input: two individuals
// output: returns true if personA has a smaller amount of fuel used than personB
bool LowerFuelSpent(const Adult& personA, const Adult& personB);

// Compare two individuals by their triptime values, used in standard sort to have array contain lowest triptime individual at start
// input: two individuals
// output: returns true if personA has a lower triptime than personB
bool LowerTripTime(const Adult& personA, const Adult& personB);

bool LowerProgress(const Adult& personA, const Adult& personB);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//TODO: below here is probably the second most confusing part of the code. I would document this in a flowchart / readme.md to draw out the particulars of how this engages with rankDistanceSort and the genetic algorithm

//Compare two adults to see if the first adult dominates the second adult
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
//  This could mean A is worse in some way from B or A and B are duplicates
bool dominationCheck(Adult& personA, Adult& personB, const cudaConstants* cConstants);

//Compare two adults by their rank
//WARNING: Using this function to sort adults will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the adults are rank 1, sorting them using this method will do nothing. 
//input: two adults
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(const Adult& personA, const Adult& personB);

//Compare two individuals by their rank and distance
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a greater distance than person B, return true
//NOTE: This is used within selectParents (specifically callSorts) and reportGeneration of optimization.cu 
//      It is used in via callSorts to put the best adults into oldAdults
//      It is used in reportGeneration to print the best performing individual 
bool rankDistanceSort(const Adult& personA, const Adult& personB);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//compares two adults to see if they have the same posDiff and speedDiff
//Input: Two adults and cConstants for the tolerances
//Output: true if they have the same posDiff and speedDiff, false if they don't
bool duplicateCheck(const Adult& personA, const Adult& personB, const cudaConstants* cConstants);

//Will take the given adult vector and find any duplicate adults
//Input:    adults - Vector of adults that will be searched through
//          cConstants - the cuda constants
//          currentAnneal - for changeing the duplicate tolerance
//Output:   Duplicate individuals within the adults array will have their error status set to duplicate
//NOTE: Ideally, this function should be used after combining old/newAdults into allAdults
//      This means the function will find duplicates within the same generation but also beween generations
//NOTE: This function does not assume the adults have been sorted
void findDuplicates (std::vector<Adult>& adults, const cudaConstants* cConstants, const double& currentAnneal);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//test version of duplicateCheck that does not use cConstants
bool duplicateCheckTest(const Adult& personA, const Adult& personB);

//Compare two adults to see if the first adult dominates the second adult
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
//used for rendezvous mission in optimization
//TEST VERSION
bool dominationCheckTest(Adult& personA, Adult& personB, int missionType);

#include "adult.cpp"

#endif