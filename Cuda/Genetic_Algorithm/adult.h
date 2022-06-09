#ifndef ADULT_H
#define ADULT_H

#include "../Genetic_Algorithm/child.h"
#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"
#include <vector>
#include <string> //used for the unitTestingRankDistanceStatusPrint method so we can get rank, distance, and status easily during unit testing

//NOTE FOR FUTURE UNDERSTANDING: the term individuals within comments is used as a general term for both adults and children.
//                               Individual most often refers to adults.
//                               A good rule of thumb: Adults will never be simulated, Children will never be sorted

// Adult is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
// Adult is inheriting startParams, errorStatus, posDiff, speedDiff, and finalPos from Child
// adult inherits the genes from child to avoid calling RK when we dont need to
struct Adult: public Child {
    // double cost;    // cost value of the individual, something that the genetic algorithm is attempting to minimize
    //TODO: probably don't need isParent, this also could be STATUS
    //bool isParent; //Indicates if this adult is already chosen as a survivor, not in use

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

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    Adult();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Input: cConstants - to access c3energy value used in param calculations
    //        childParameters - struct returned by generateNewIndividual()
    // Output: this adult's startParams.y0 is set to the initial position and velocity of the spacecraft
    Adult(rkParameters<double> & childParameters, const cudaConstants* cConstants): Child(childParameters, cConstants), rank(INT_MAX), distance(-1){} 

    // Constructor to create an adult from a child
    // Input: A child - this has cConstants, a set of rkParameters, and elements in it already, as well as the speedDiff and posDiff
    // Output: a child is turned into an adult
    Adult(const Child& c): Child(c), rank(INT_MAX), distance(-1){}

    // Constructor that will create an adult from another adult
    // Input: an adult; the variables for the adult passed in will be copied over to the new adult
    // Output: a new adult, which is a copy of the passed-in adult
    //      NOTE: this is needed for the function of shuffling the oldAdult vector within the newGeneration function of ga_crossover.cpp.
    //            This is needed to create the temp Adult needed to do the swap
    Adult(const Adult& a): Child(a), rank(a.rank), distance(a.distance){}

    //TODO: Consider deleting this -> it is weird sorting from best to worst and calling best less than worst (unless we want to be very specific we mean LESS WRONG or something?)
    //Compare two adults by their rank and distance
    //input: another adult
    //output: if this adult's rank is lower than the other adult's rank, return true
    //        if this adult and the other adult have the same rank and this adult has a greater distance than the other adult, return true
    //TODO: What about status? We should probably document that here
    //Sorts the whole pool from lowest to highest rank. Adults of the same rank are sorted from highest to lowest distance
    bool operator<(const Adult &other);

    //TODO: Comments on where these ^^^^ constructors should be called (i.e. at the end of newGeneration, after callRK)

    //#ifdef UNITTEST //this should be defined in unit testing

    // Constructor used ONLY FOR UNIT TESTING!!!
    // This constructor uses the default rkParameters, elements, posDiff and speedDiff, the only things that are defined are the rank and distance because that is all rankDistance sort cares about
    // Input: A rank (r), a distance (d), and a STATUS (s) -> status defaults to VALID if undefined
    // Output: Constructs an Adult using the unit testing Child constructor that only contains default and sets its rank and distance so it can be sorted using rankDistanceSort
    // ONLY USE FOR UNIT TESTING!!! SHOULD NEVER BE CALLED ELSEWHERE 
    Adult(int r, double d, int s = VALID): Child(s), rank(r), distance(d){}
   
    // A function used exclusively for unit testing that prints that rank and distance of an adult
    // Output: takes an adult's rank, distance, and status and creates a string holding this information
    // ONLY USE FOR UNIT TESTING!!! Note: This was purposefully given a long name to discourage people from using it
    const std::string unitTestingRankDistanceStatusPrint();

    //getters used exclusively in unit testing
    int getRank();
    double getDistance();
    //#endif //unit testing endif

};

//TODO: document style: where are these being used? They are pretty major, we should help clarify

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

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two adults to see if the first adult dominates the second adult
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
//used for rendezvous mission in optimization
bool dominationCheck(Adult& personA, Adult& personB, const cudaConstants* cConstants);

//Compare two adults by their rank
//WARNING: Using this function to sort adults will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the adults are rank 1, sorting them using this method will do nothing. 
//input: two adults
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(const Adult& personA, const Adult& personB);

//Compare two individuals by their rank and distance
//input: two individuals
//TODO: consider status (this goes for the above ones)
//TODO: Where is this called ? Used?
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a greater distance than person B, return true
bool rankDistanceSort(const Adult& personA, const Adult& personB);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two adults to see if the first adult dominates the second adult
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
//used for rendezvous mission in optimization
//TEST VERSION
bool dominationCheckTest(Adult& personA, Adult& personB);

#include "adult.cpp"

#endif