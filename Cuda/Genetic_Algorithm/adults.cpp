#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI
#include "adults.h"

// Default constructor
Adult::Adult(){
    dominatedByCount = 0;
    rank = 0; //might change later?
    distance = -1;
    isParent = false;
}


//Compare two adults by their rank and distance
//input: another adult
//output: if this adult's rank is lower than the other adult's rank, return true
//        if this adult and the other adult have the same rank and this adult has a greater distance than the other adult, return true
//Sorts the whole pool from lowest to highest rank. Adults of the same rank are sorted from highest to lowest distance
bool Adult::operator<(const Adult &other) {
    //TODO: is this a good system to check validity first?
    if (status != VALID){
        return false;
    }
    if(rank < other.rank){
        return true;
    }
    else if (rank == other.rank && distance > other.distance){
        return true;
    }
    else {
        return false;
    }
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Compare two individuals by their positional difference values, used in standard sort to have array contain lowest posDiff individual at start
// input: two individuals
// output: returns true if personB has a higher positional difference than personA
bool LowerPosDiff(Adult& personA, Adult& personB) {
    if (personA.posDiff < personB.posDiff) {
        return true;
    }
    else {
        return false;
    }

}

// Compare two individuals by their velocity difference values, used in standard sort to have array contain highest speedDiff individual at start
// input: two individuals
// output: returns true if personB has a lower velocity difference than personA
bool HigherSpeedDiff(const Adult& personA, const Adult& personB) {
    if (personA.speedDiff > personB.speedDiff) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two individuals by their velocity difference values, used in standard sort to have array contain lowest speedDiff individual at start
// input: two individuals
// output: returns true if personA has a lower velocity difference than personB
bool LowerSpeedDiff(const Adult& personA, const Adult& personB) {
    if (personA.speedDiff < personB.speedDiff) {
        return true;
    }
    else {
        return false;
    }

}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals to see if the first individual dominates the second individual
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
bool dominates(Adult& personA, Adult& personB, const cudaConstants* cConstants) {

    //TODO: Might want to consider modifying tolerances 
    
    //Is true if A is at least equally as good as B for all objectives
    bool AisEqual = false;
    //Is true if A is better than B for at least one objective
    bool AisBetter = false;
    //tolerances used to determine the range of values considered equal
    double posTolerance = cConstants->posDominationTolerance;
    double speedTolerance = cConstants->speedDominationTolerance;
    //true if A posDiff "equals" B posDiff
    bool APosEqualsB = false;
    //true if A speedDiff "equals" B speedDiff
    bool ASpeedEqualsB = false;

    //True is A posdiff is equal to B posDiff +- posTolerance
    if ((personA.posDiff < personB.posDiff + posTolerance) && (personA.posDiff > personB.posDiff - posTolerance)){
        APosEqualsB = true;
    }
    //True is A speeddiff is equal to B speedDiff +- speedTolerance
    if ((personA.speedDiff < personB.speedDiff + speedTolerance) && (personA.speedDiff > personB.speedDiff - speedTolerance)){
        ASpeedEqualsB = true;
    }
    //If A.posDiff is approximately/better than B.posDiff, and A.speedDiff is approximately/better than B.speedDiff, then A is equal to B.
    if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff < personB.speedDiff || ASpeedEqualsB)) {
        AisEqual = true;
    }
    //If A has a better posDiff or speedDiff than B, then A is better than B
    if (personA.posDiff < personB.posDiff || personA.speedDiff < personB.speedDiff){
        AisBetter = true;
    }

    //A Dominates B
    //A only dominates B if:
        //A.posDiff - posTolerance < B.posDiff <= A.posDiff, and A.speedDiff < B.speedDiff
        //A.posDiff < B.posDiff, and A.speedDiff - speedTolerance < B.speedDiff <= A.speedDiff
        //A.posDiff < B.posDiff and A.speedDiff < B.speedDiff (better in every way)
    if (AisEqual && AisBetter){
        return true;
    }
    //A and B are codominant if:
        //A.posDiff < B.posDiff & A.speedDiff > B.speedDiff
        //A.posDiff > B.posDiff & A.speedDiff < B.speedDiff
        //A.posDiff = B.posDiff & A.speedDiff = B.speedDiff (Unlikely since they're doubles, unless A & B are clones or genetically similar)
    //B dominates A if any of the conditions for "A dominates B", when reversed, are true.
    else {
        return false;
    }
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals by their rank
//WARNING: Using this function to sort individuals will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the individuals are rank 1, sorting them using this method will do nothing. 
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(const Adult& personA, const Adult& personB){
    //TODO: is this a good system to check validity first?
    if(personA.status != VALID){ //if personA has nan values or other errors they are set as worse than other adults (even other adults with errors)
        return false;
    }
    if (personA.rank < personB.rank) {
        return true;
    }
    else {
        return false;
    }
}


//Compare two individuals by their rank and distance
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a greater distance than person B, return true
//Sorts the whole pool from lowest to highest rank. Individuals of the same rank are sorted from highest to lowest distance
bool rankDistanceSort(const Adult& personA, const Adult& personB) {
    //TODO: is this a good system to check validity first?
    if (personA.status != VALID){
        return false;
    }
    if(personA.rank < personB.rank){
        return true;
    }
    else if (personA.rank == personB.rank && personA.distance > personB.distance){
        return true;
    }
    else {
        return false;
    }
    // OLD VERSION....
    // if(personA.rank < personB.rank){
    //     return true;
    // }
    // else if (personA.rank == personB.rank && personA.distance > personB.distance){
    //     return true;
    // }
    // else {
    //     return false;
    // }
    

}



