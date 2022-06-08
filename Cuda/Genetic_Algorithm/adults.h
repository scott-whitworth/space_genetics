#ifndef ADULTS_H
#define ADULTS_H

#include "../Genetic_Algorithm/children.h"
#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"
#include <vector>

//TODO: Clarify / change 'individual' term vs. 'Child' and 'Adult'

//TODO: documentation to be VERY clear why we are inheriting 
//      gene-perspective would be helpful
// Individual is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
struct Adult: public Child {
   // double cost;    // cost value of the individual, something that the genetic algorithm is attempting to minimize
   //TODO: probably don't need isParent, this also could be STATUS
    bool isParent; //Indicates if this individual is already chosen as a survivor

    int dominatedByCount; //keeps track of how many individuals this individual has been dominated by
    //TODO: more documentation (like this is used by the dominated sort...)
    int rank; //how good this individual is compared to other individuals in it's generation. The lower the rank, the better.

    std::vector<int> dominates; //The pool indexes of individuals that have been dominated by this individual 
    //TODO: See comment in optimiaztion.cu, giveRank() about moving dominates and dominatedByCount
    double distance; //how different this individual is from the individuals most similar to it.

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    Adult();

    //TODO: Rank 0 might not be the best value to initialize rank to - consider whether or not this should be changed
    //      This is probably just assigning max int (google search: c++ MAX_INT value)

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Input: cConstants - to access c3energy value used in getCost()
    //        newInd - struct returned by generateNewIndividual()
    // Output: this adult's startParams.y0 is set to the initial position and velocity of the spacecraft
    Adult(rkParameters<double> & newInd, const cudaConstants* cConstants): Child(newInd, cConstants), dominatedByCount(0), rank(0), distance(-1), isParent(false){}

    // Constructor to create an adult from a child
    // Input: A child - this has cConstants, a set of rkParameters, and elements in it already, as well as the speedDiff and posDiff
    // Output: a child is turned into an adult
    Adult(const Child& c): Child(c), dominatedByCount(0), rank(0), distance(-1), isParent(false){}

    //TODO: Comments on where these ^^^^ constructors should be called (i.e. at the end of newGeneration, after callRK)

    
    //Compare two adults by their rank and distance
    //input: another adult
    //output: if this adult's rank is lower than the other adult's rank, return true
    //        if this adult and the other adult have the same rank and this adult has a greater distance than the other adult, return true
    //TODO: What about status? We should probably document that here
    //Sorts the whole pool from lowest to highest rank. Adults of the same rank are sorted from highest to lowest distance
    bool operator<(const Adult &other);

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
//Compare two adults to see if the first individual dominates the second adult
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
bool dominates(Adult& personA, Adult& personB, const cudaConstants* cConstants);

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

#include "adults.cpp"

#endif