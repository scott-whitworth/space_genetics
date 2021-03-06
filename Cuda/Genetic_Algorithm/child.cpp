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
    testCount = 0;
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
    fuelSpent = other.fuelSpent;
    funcStatus = other.funcStatus;
    errorStatus = other.errorStatus;
    birthday = other.birthday;
    avgParentProgress = other.avgParentProgress;
    progress = other.progress;
    testCount = other.testCount;
}

//Getter for a parameter dependent on the objective that is passed in
__host__ double Child::getParameters (const objective & requestObjective) const {
    //if/esle tree will go find the parameter goal of the request objective and return the associated value
    if (requestObjective.goal == MIN_POS_DIFF) {
        return posDiff;
    }
    else if (requestObjective.goal == MIN_SPEED_DIFF || requestObjective.goal == MAX_SPEED_DIFF) {
        return speedDiff;
    }
    else if (requestObjective.goal == MIN_FUEL_SPENT) {
        return fuelSpent;
    }
    else if (requestObjective.goal == MIN_TRIP_TIME) {
        return startParams.tripTime; 
    }
    else if (requestObjective.goal == MIN_ORBIT_POS_DIFF){
        return posDiff;
    }
    else if (requestObjective.goal == MIN_ORBIT_SPEED_DIFF){
        return speedDiff;
    }
    else {
        //Inducates error
        std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        return -1;
    }
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_target, theta_fin_target, and z_fin_target
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Child::getPosDiff(const cudaConstants* cConstants) {
    
    
    //Check to see if the posDiff should be calculated or set to the bad value
    //Will mean only valid children will be likely to be considered for future generations
    if (errorStatus == VALID || errorStatus == DUPLICATE) {
        //Check to see if getting into an orbit is the goal
        if (cConstants->orbitalRadius != NOT_APPLICABLE) {
            //Calculate the orbital posDiff
            //Calculates the distance between the target and the spacecraft
            double targetToCraftDist = (sqrt(pow(cConstants->r_fin_target, 2) + pow(finalPos.r, 2) + pow(finalPos.z - cConstants->z_fin_target, 2) - (2*finalPos.r*cConstants->r_fin_target*cos(cConstants->theta_fin_target-finalPos.theta))));
            //The distance between the target and the spacecraft should be the orbital radius in order to converge
            posDiff = abs(targetToCraftDist - cConstants->orbitalRadius);
        }
        else {
            //Calculate the normal posDiff
            //posDiff = sqrt(delta(r)^2 + delta(theta)^2 - r*fmod(theta, 2pi) + delta(z)^2 ) -> magnitude of delta(position)
            posDiff = sqrt(pow(cConstants->r_fin_target - finalPos.r, 2) + pow( (cConstants->r_fin_target * cConstants->theta_fin_target) - (finalPos.r * fmod(finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_target - finalPos.z, 2));
        }
    }

    return posDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_target, vtheta_fin_target, and vz_fin_target
// Output: Assigns and returns this child's speedDiff value
__host__ __device__ double Child::getSpeedDiff(const cudaConstants* cConstants) {
    //Check to see if the speedDiff should be calculated or set to the bad value
    //Will mean only valid children will be likely to be considered for future generations
    if (errorStatus == VALID || errorStatus == DUPLICATE) {
        //Check to see if the goal is an orbit vs another mission type
        if (cConstants->orbitalSpeed != NOT_APPLICABLE){
            //Goal is an orbit
            // The radial velocity of the spacecraft compared with that of the target should be constant to remain in a circular orbit
            //      Thus the pow(cConstants->vr_fin_target - finalPos.vr, 2) should be 0
            // Additionally, to remain within a circular orbit, the magnitude of the sqrt(v_theta ^2 + v_z^2) must be the orbital velocity
            //      So the (sqrt(v_theta ^2 + v_z^2) - v_orbit)^2 term shoudl also go to 0
            speedDiff = sqrt(pow(cConstants->vr_fin_target - finalPos.vr, 2) + pow(sqrt(pow(cConstants->vtheta_fin_target - finalPos.vtheta, 2) + pow(cConstants->vz_fin_target - finalPos.vz, 2)) - cConstants->orbitalSpeed, 2));
        }
        else{
            //Goal is impact or rendezvous
            speedDiff = sqrt(pow(cConstants->vr_fin_target - finalPos.vr, 2) + pow(cConstants->vtheta_fin_target - finalPos.vtheta, 2) + pow(cConstants->vz_fin_target - finalPos.vz, 2)); 
        }
    }

    return speedDiff;
}

__host__ void Child::getProgress(const cudaConstants* cConstants) {

    //Holds the progress of the child before it is actually assigned to the child
    //  The relative costs of each objective will be added to it
    //  It will then be divided by the number of objectives and assigned to the child
    double calcProgress = 0; 

    //Iterate through the objectives
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
        //Check to see if the goal for this objective is to minimize or maximize the parameter 
        //  Necessary because the progress values are calculated differently depending on the direction 
        //  See the objective header for details on how objective direction is determined
        if (cConstants->missionObjectives[i].goal < 0) {//Minimization
            
            //See if the child has met the the convergence threshold for this parameter
            if (getParameters(cConstants->missionObjectives[i]) < cConstants->missionObjectives[i].convergenceThreshold) {
                //Add one to the progress to signify that the parameter has met the goal
                calcProgress += 1; 
            }
            //The child hasn't met the parameter goal
            else {
                //Add the progress for this parameter to the goal
                //For minimization, the progress is the parameter divided by the threshold
                calcProgress += (getParameters(cConstants->missionObjectives[i])/cConstants->missionObjectives[i].convergenceThreshold); 
            }
        }
        //Maximization is very similar minimization, but the signs are flipped and an inverse fraction is used
        else if (cConstants->missionObjectives[i].goal > 0) {//Maximization
            
            //See if the child has met the the convergence threshold for this parameter
            if (getParameters(cConstants->missionObjectives[i]) > cConstants->missionObjectives[i].convergenceThreshold) {
                //Add one to the progress to signify that the parameter has met the goal
                calcProgress += 1; 
            }
            //The child hasn't met the parameter goal
            else {
                //Add the progress for this parameter to the goal
                //For maximization, the progress is the threshold divided by the parameter
                calcProgress += (cConstants->missionObjectives[i].convergenceThreshold/getParameters(cConstants->missionObjectives[i])); 
            }
        }
        //No mission type was identified 
        else {
            std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        }
    }

    //The total cost has been calculated
    //It needs to be divided by the number of objectives to find the weighted average progress for each objective
    calcProgress = cConstants->missionObjectives.size()/calcProgress;

    //Assign the weighted progress to the child
    progress = calcProgress; 
}

//TODO: Long term we need to consider the velocity vector, not just the magnitude
//      this would be where we want to introduce those calculations
//      Three parameters:
//              position diff
//              vector direction (unit vector) diff (this is what we will need to implement)
//              magnitude of velocity diff (this is what we have)
