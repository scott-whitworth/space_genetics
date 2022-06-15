#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
//this is never called, and if it is we will know because of the func status
Child::Child() {
    //TODO: all we do is set status to DEFAULT_CHILD (which is essentially an error if we ever try to process)
    //TODO: get rid of posDiff/speedDiff
    //posDiff = 1.0;
    //speedDiff = 0.0; //TODO: This is ok for impact, but an issue for soft landing

    funcStatus = DEFAULT_CHILD;//not ready to be an adult

    errorStatus = NOT_RUN; //not run through callRk
}

// Set the initial position of the spacecraft according to the newly generated parameters
// Input: cConstants - to access c3energy value used in getCost()
//        childParameters - struct returned by generateNewIndividual()
// Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
Child::Child(rkParameters<double> & childParameters, const cudaConstants* cConstants) {

    startParams = childParameters;
    elements<double> earth = launchCon->getCondition(startParams.tripTime); //get Earth's position and velocity at launch

    startParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(startParams.alpha),
        earth.theta+asin(sin(M_PI-startParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(startParams.zeta)*sin(startParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(startParams.zeta)*cos(startParams.beta)*cConstants->v_escape,
        earth.vz+sin(startParams.zeta)*cConstants->v_escape);

    funcStatus = FUNCTIONAL_CHILD;//ready to be an adult
    errorStatus = NOT_RUN; //not run through callRK yet
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
    funcStatus = other.funcStatus;
    errorStatus = other.errorStatus;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Child::getPosDiff(const cudaConstants* cConstants) {
    //Check to see if the posDiff should be calculated or set to the bad value
    //Will mean only valid children will be likely to be considered for future generations
    if (errorStatus == VALID) {
        //posDiff = sqrt(delta(r)^2 + delta(theta)^2 - r*fmod(theta, 2pi) + delta(z)^2 ) -> magnitude of delta(position)
        posDiff = sqrt(pow(cConstants->r_fin_ast - finalPos.r, 2) + pow( (cConstants->r_fin_ast * cConstants->theta_fin_ast) - (finalPos.r * fmod(finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_ast - finalPos.z, 2));
    }
    else {
        posDiff = BAD_POSDIFF;
    }

    return posDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
// Output: Assigns and returns this child's speedDiff value
__host__ __device__ double Child::getSpeedDiff(const cudaConstants* cConstants) {
    //Check to see if the speedDiff should be calculated or set to the bad value
    //Will mean only valid children will be likely to be considered for future generations
    if (errorStatus == VALID) {
        //TODO: Check if this is accurate:
        //speedDiff = sqrt(delta(vr)^2 + delta(vtheta)^2 + delta(vz)^2 ) -> magnitude of delta(position)
        //we want to include: v_ast = sqrt(pow(cConstants->vr_fin_ast, 2) + pow(cConstants->vtheta_fin_ast, 2) + pow(cConstants->vz_fin_ast, 2))?
        //                    v_pos = sqrt(pow(finalPos.vr, 2) + pow(finalPos.vtheta, 2) + pow(finalPos.vz, 2))?
        //                    vDiff = abs(v_ast - v_pos)?
        speedDiff = sqrt(pow(cConstants->vr_fin_ast - finalPos.vr, 2) + pow(cConstants->vtheta_fin_ast - finalPos.vtheta, 2) + pow(cConstants->vz_fin_ast - finalPos.vz, 2)); 
    }
    else {
        //The bad speedDiff is dependent on the type of mission; a high speedDiff is bad for rendezvous missions and a low speedDiff is bad for impact missions
        if (cConstants -> missionType == Impact)
        {
            speedDiff = BAD_IMPACT_SPEEDDIFF;
        }
        else {
            speedDiff = BAD_RENDEV_SPEEDDIFF;
        }
        
    }
    return speedDiff;
}

//TODO: Long term we need to consider the velocity vector, not just the magnitude
//      this would be where we want to introduce those calculations
//      Three parameters:
//              position diff
//              vector direction (unit vector) diff (this is what we will need to implement)
//              magnitude of velocity diff (this is what we have)
