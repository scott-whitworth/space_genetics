#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
Individual::Individual() {
    this->posDiff = 1.0;
    this->speedDiff = 0.0; //TODO: This is ok for impact, but an issue for soft landing
    this->cost = 10;
    this->dominatedCount = 0;
    this->rank = 0; //might change later?
    this->distance = -1;
    this->isParent = false;
}

// Set the initial position of the spacecraft according to the newly generated parameters
// Input: cConstants - to access c3energy value used in getCost()
//        newInd - struct returned by generateNewIndividual()
// Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
Individual::Individual(rkParameters<double> & newInd, const cudaConstants* cConstants) {

    this->startParams = newInd;
    elements<double> earth = launchCon->getCondition(this->startParams.tripTime); //get Earth's position and velocity at launch

    this->startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(this->startParams.alpha),
        earth.theta+asin(sin(M_PI-this->startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(this->startParams.zeta)*sin(this->startParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(this->startParams.zeta)*cos(this->startParams.beta)*cConstants->v_escape,
        earth.vz+sin(this->startParams.zeta)*cConstants->v_escape);
    
    this->cost = 10;
    this->dominatedCount = 0;
    this->rank = 0; //might change later?
    this->distance = -1;
    this->isParent = false;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Individual::getPosDiff(const cudaConstants* cConstants) {
   this->posDiff = sqrt(pow(cConstants->r_fin_ast - this->finalPos.r, 2) + pow( (cConstants->r_fin_ast * cConstants->theta_fin_ast) - (this->finalPos.r * fmod(this->finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_ast - this->finalPos.z, 2));
   return this->posDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
// Output: Assigns and returns this individual's speedDiff value
__host__ __device__ double Individual::getSpeedDiff(const cudaConstants* cConstants) {
    this->speedDiff = sqrt(pow(cConstants->vr_fin_ast - this->finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - this->finalPos.vtheta, 2) + pow(cConstants->vz_fin_ast - this->finalPos.vz, 2)); 
    return this->speedDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a cost value to quantitatively evaluate this Individual
// Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
// Output: Assigns and returns this individuals cost value
// TODO: This is not final, could change / be optimized

// __host__ __device__ double getCost(const cudaConstants* cConstants) {
//     this->cost = this->posDiff;
//     return this->cost;
// }

__host__ __device__ double Individual::getCost_Hard(const cudaConstants* cConstants) {
    this->cost = this->posDiff;//This is in AU
    return this->cost;
}

__host__ __device__ double Individual::getCost_Soft(const cudaConstants* cConstants) {
    //this->cost = sqrt(pow((this->posDiff), 2) + pow((this->speedDiff*cConstants->timeRes), 2));//This is in AU when speedDiff is multiplied by seconds
    double posProximity = this->posDiff - cConstants->pos_threshold;
    double speedProximity = this->speedDiff - cConstants->speed_threshold;

    //If it's already passed the tolerance in one objective, don't let it become negative
    if (posProximity < 0) {
        posProximity = 0;
    }
    if (speedProximity < 0) {
        speedProximity = 0;
    }
    
    this->cost = posProximity + speedProximity; //TODO: This may need to change
    
    return this->cost;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Individual::operator>(Individual &other) {
    //TODO:: This should call getCost
    if (this->cost > other.cost) {
        return true;
    }
    else {
        return false;
    }
}


bool Individual::operator<(Individual &other) {
    if (this->cost < other.cost) {
        return true;
    }
    else {
        return false;
    }
}

bool Individual::operator==(Individual &other) {
    if (this->cost == other.cost) {
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

//TODO:: These should be const references
bool LowerPosDiff(Individual& personA, Individual& personB) {
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
bool HigherSpeedDiff(Individual& personA, Individual& personB) {
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
bool LowerSpeedDiff(Individual& personA, Individual& personB) {
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
bool dominates(Individual& personA, Individual& personB, const cudaConstants* cConstants) {

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
bool rankSort(Individual& personA, Individual& personB){

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
bool rankDistanceSort(Individual& personA, Individual& personB) {

    if(personA.rank < personB.rank){
        return true;
    }
    else if (personA.rank == personB.rank && personA.distance > personB.distance){
        return true;
    }
    else {
        return false;
    }
}