#include <math.h>
#include "children.h"
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
Child::Child() {
    posDiff = 1.0;
    speedDiff = 0.0; //TODO: This is ok for impact, but an issue for soft landing
}

// Set the initial position of the spacecraft according to the newly generated parameters
// Input: cConstants - to access c3energy value used in getCost()
//        newInd - struct returned by generateNewIndividual()
// Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
Child::Child(rkParameters<double> & newChild, const cudaConstants* cConstants) {

    startParams = newChild;
    elements<double> earth = launchCon->getCondition(startParams.tripTime); //get Earth's position and velocity at launch

    startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(startParams.alpha),
        earth.theta+asin(sin(M_PI-startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(startParams.zeta)*sin(startParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(startParams.zeta)*cos(startParams.beta)*cConstants->v_escape,
        earth.vz+sin(startParams.zeta)*cConstants->v_escape);
}

// Copy constructor
// Sets this child's parameters, elements, posDiff, and speedDiff to another child's values for these quantities
// Input: Another child
// Output: this child is a copy of the other child that was passed in
Child:: Child(const Child& other){
    startParams = other.startParams;
    finalPos = other.finalPos;
    posDiff = other.posDiff; 
    speedDiff = other.speedDiff;
    status = other.status;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Child::getPosDiff(const cudaConstants* cConstants) {
   posDiff = sqrt(pow(cConstants->r_fin_ast - finalPos.r, 2) + pow( (cConstants->r_fin_ast * cConstants->theta_fin_ast) - (finalPos.r * fmod(finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_ast - finalPos.z, 2));
   return posDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
// Output: Assigns and returns this child's speedDiff value
__host__ __device__ double Child::getSpeedDiff(const cudaConstants* cConstants) {
    speedDiff = sqrt(pow(cConstants->vr_fin_ast - finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - finalPos.vtheta, 2) + pow(cConstants->vz_fin_ast - finalPos.vz, 2)); 
    return speedDiff;
}
