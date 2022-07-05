#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
//this is never called, and if it is we will know because of the func status
Child::Child() {

    //posDiff = 1.0;
    //speedDiff = 0.0; //This is ok for impact, but an issue for soft landing

    funcStatus = DEFAULT_CHILD;//not ready to be an adult

    errorStatus = NOT_RUN; //not run through callRk
}

// Set the initial position of the spacecraft according to the newly generated parameters
// Input:  childParameters - the starting paramters that are generated randomly, from mutation, or from crossover
//         cConstants - to access the v_escape value for calculations of the startParams
//         genCreated - the generation this child is created, used for birthday/age data
//         calcAvgParentProgress - the creating parent's average progress. Will be set to 0 for randomly generated children
// Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
Child::Child(rkParameters<double> & childParameters, const cudaConstants* cConstants, int genCreated, double calcAvgParentProgress) {

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
    birthday = genCreated; //Set the child's birthday to the current generation
    avgParentProgress = calcAvgParentProgress; //The avg progress of the creating parents, if any (0 for randomly generated children)
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
    birthday = other.birthday;
    avgParentProgress = other.avgParentProgress;
    progress = other.progress;

}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Child::getPosDiff(const cudaConstants* cConstants) {
    //Check to see if the posDiff should be calculated or set to the bad value
    //Will mean only valid children will be likely to be considered for future generations
    if (errorStatus == VALID || errorStatus == DUPLICATE) {
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
    if (errorStatus == VALID || errorStatus == DUPLICATE) {
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

__host__ __device__ void Child::getProgress(const cudaConstants* cConstants) {

    //Variables will be used to calculate progress
    //      They both will specifically be used to calculate the adult's diff / goal diff
    //      This prevents a excellent speed diff causing the adult's diff / goal diff to be less than 1
    //          If this was the case, there is the potential that the progress goal could be met even if the position diff hadn't hit it's goal
    double rel_posDiff, rel_speedDiff;
    rel_posDiff = rel_speedDiff = 0;

    //How to calculate the progress will depend on the mission type
    //In general, progress will be the number of objectives / the combined rel diffs of the objectives
    //      The rel diffs are themselves calculated by dividing the current diff by the threshold of that objective, which calculates how close the child is to that goal
    //      To make sure an individual beating one goal by a lot but not meeting the other doesn't have a progress of 1, if the diffs are meeting the goal, the rel diffs are set to 1 to signify 100% progress
    //This will mean a progress of 1 signifies that the individual has hit the mission goals

    //For impacts, the progress only depends on posDiff
    //  NOTE: While the RD sorting favors high speed for impacts, convergence ultimately comes down to posDiff
    if (cConstants -> missionType == Impact) {
        //The goal is position threshold and the current status is the best individual's posDiff

        //Check to see if the rel_posDiff would be less than the threshold
        if (posDiff > cConstants->pos_threshold) {
            //If not, calculate how close the diff is to the goal
            rel_posDiff = posDiff / cConstants->pos_threshold; 
        }
        else {
            //If so, set the rel_posDiff to 1, signifying that the goal has been met
            rel_posDiff = 1.0;
        }

        //The progress is measured by the number of objectives (1) divided by the rel diff
        //So, a simulation meeting the requirements has a progress value of 1 since it would be 1/1
        progress = Impact/rel_posDiff;
    }
    //For rendezvous, the progress depends on both posDiff and speedDiff
    else {
        //Similarly to impact, calculate how far the best adult's speed and position diffs are away from the goal and subtract by the number of mission goals
        
        //First, same as above either set rel_posDiff to how close the best adult's posDiff is to the goal or to 1, depending on if the posDiff is beating the goal
        if (posDiff > cConstants->pos_threshold) {
            rel_posDiff = posDiff / cConstants->pos_threshold;
        }
        else {
            rel_posDiff = 1.0;
        }
        //Repeat the same process in calculating rel_posDiff for rel_speedDiff
        if (speedDiff > cConstants->speed_threshold) {
            rel_speedDiff = speedDiff / cConstants->speed_threshold;
        }
        else { 
            rel_speedDiff = 1.0;
        }

        //The progress calculates as the number of objectives for the mission (2) divided by the sum of both rel diffs
        //This means a flight that meets all goals will have a progress of 1, since the rel pos/speed diffs would be set to 1, so it would be 2/2
        progress = Rendezvous/(rel_posDiff + rel_speedDiff);
    }
    
    //End the function
    return; 
}

//TODO: Long term we need to consider the velocity vector, not just the magnitude
//      this would be where we want to introduce those calculations
//      Three parameters:
//              position diff
//              vector direction (unit vector) diff (this is what we will need to implement)
//              magnitude of velocity diff (this is what we have)
