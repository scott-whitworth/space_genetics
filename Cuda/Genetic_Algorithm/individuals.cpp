#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
Individual::Individual() {
    this->posDiff = 1.0;
    this->speedDiff = 0.0;
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
}

// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Individual::getPosDiff(const cudaConstants* cConstants) {
    this->posDiff = sqrt(pow(cConstants->r_fin_ast - this->finalPos.r, 2) + pow( (cConstants->r_fin_ast * cConstants->theta_fin_ast) - (this->finalPos.r * fmod(this->finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_ast - this->finalPos.z, 2));
//    this->posDiff = calcPosDiff(this->finalPos, cConstants);
   return this->posDiff;
}

//Utility to calculate the position difference between a position velocity set and the asteroid final position
//__host__ __device__ double calcPosDiff(const elements<double>& currentState, const cudaConstants* cConstants){
//    return  sqrt(pow(cConstants->r_fin_ast - currentState.r, 2) + pow( (cConstants->r_fin_ast * cConstants->theta_fin_ast) - (currentState.r * fmod(currentState.theta, 2 * M_PI) ) , 2) + pow(cConstants->z_fin_ast - currentState.z, 2));
//}

// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
// Output: Assigns and returns this individual's speedDiff value
__host__ __device__ double Individual::getSpeedDiff(const cudaConstants* cConstants) {
//    this->speedDiff = sqrt(pow(cConstants->vr_fin_ast - this->finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - this->finalPos.vtheta, 2) + pow(cConstants->vz_fin_ast - this->finalPos.vz, 2));
    this->speedDiff = sqrt(pow(cConstants->vr_fin_ast - this->finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - this->finalPos.vtheta, 2) + pow(cConstants->vz_fin_ast - this->finalPos.vz, 2)); 
    return this->speedDiff;
}

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
    this->cost = sqrt(pow(this->posDiff, 2) + pow((cConstants->timeRes*(this->speedDiff)), 2));//This is in AU
    return this->cost;
}

bool Individual::operator>(Individual &other) {
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

// Compare two individuals by their positional difference values, used in standard sort to have array contain lowest posDiff individual at start
// input: two individuals
// output: returns true if personB has a higher positional difference than personA
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

