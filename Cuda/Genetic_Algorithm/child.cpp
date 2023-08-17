#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI


// Default constructor
//this is never called, and if it is we will know because of the func status
Child::Child() {

    errorStatus = NOT_RUN; //not run through callRk

    simStatus = INITIAL_SIM; //has not been simulated
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

    errorStatus = NOT_RUN; //not run through callRK yet
    simStatus = INITIAL_SIM; //Has not been simulated yet

    birthday = genCreated; //Set the child's birthday to the current generation

    avgParentProgress = calcAvgParentProgress; //The avg progress of the creating parents, if any (0 for randomly generated children)
    normalizedObj = std::vector<double>(cConstants->missionObjectives.size(), 0); //Resize the objective progress vector to match the number of objectives and set the default progress to be 0 (bad value)

    stepCount = 0; //no calculations done yet, default is zero
    simStartTime = 0; //Inititially start the simulation at the start of the trip time
    simNum = 0; //Has not been simulated yet

    minMarsDist = 100; //Arbitrarily high initial min mars distance
    orbithChange = 1e-14; //No assist initially, so no angular momentum change initally
}

// Copy constructor (Needed for sorts, so make sure this is up to date with the child's variables!)
// Sets this child's parameters, elements, posDiff, and speedDiff to another child's values for these quantities
// Input: Another child
// Output: this child is a copy of the other child that was passed in
Child:: Child(const Child& other){
    startParams = other.startParams;
    finalPos = other.finalPos;
    posDiff = other.posDiff; 
    speedDiff = other.speedDiff;
    fuelSpent = other.fuelSpent;
    orbitPosDiff = other.orbitPosDiff;
    orbitSpeedDiff = other.orbitSpeedDiff; 
    errorStatus = other.errorStatus;
    simStatus = other.simStatus;
    birthday = other.birthday;
    avgParentProgress = other.avgParentProgress;
    progress = other.progress;
    normalizedObj = other.normalizedObj;
    stepCount = other.stepCount;
    minMarsDist = other.minMarsDist;
    orbithChange = other.orbithChange;
    simNum = other.simNum;
}

//Getter for a parameter dependent on the objective that is passed in
__host__ double Child::getParameters (const objective & requestObjective) const {
    //if/esle tree will go find the parameter goal of the request objective and return the associated value
    if (requestObjective.goal == MIN_POS_DIFF) {
        return posDiff;
    }
    else if (requestObjective.goal == MIN_SPEED_DIFF) {
        return speedDiff;
    }
    else if (requestObjective.goal == MAX_SPEED_DIFF) {
        return 1/speedDiff;
    }
    else if (requestObjective.goal == MIN_FUEL_SPENT) {
        return fuelSpent;
    }
    else if (requestObjective.goal == MIN_TRIP_TIME) {
        return startParams.tripTime; 
    }
    else if (requestObjective.goal == MIN_ORBIT_POS_DIFF){
        return orbitPosDiff;
    }
    else if (requestObjective.goal == MIN_ORBIT_SPEED_DIFF){
        return orbitSpeedDiff;
    }
    else if (requestObjective.goal == MIN_MARS_DIST){
        return minMarsDist;
    }
    else if (requestObjective.goal == MAX_ORBIT_ASST){
        return 1/orbithChange;
    }
    else {
        //Indicates error
        std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        return -1;
    }
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Calculates a posDiff value
// Input: cConstants in accessing properties such as r_fin_target, theta_fin_target, and z_fin_target
// Output: Assigns and returns this individual's posDiff value
__host__ __device__ double Child::getPosDiff(const cudaConstants* cConstants) {
    //Calculate the normal posDiff using the sun-relative pos diff between the spacecraft and the target 
    //posDiff = sqrt(delta(r)^2 + delta(theta)^2 - r*fmod(theta, 2pi) + delta(z)^2 ) -> magnitude of delta(position)
    posDiff = sqrt(pow(cConstants->r_fin_target - finalPos.r, 2) + pow( (cConstants->r_fin_target * cConstants->theta_fin_target) - (finalPos.r * fmod(finalPos.theta, 2 * M_PI)), 2) + pow(cConstants->z_fin_target - finalPos.z, 2));

    return posDiff;
}

// Calculates a speedDiff value
// Input: cConstants in accessing properties such as vr_fin_target, vtheta_fin_target, and vz_fin_target
// Output: Assigns and returns this child's speedDiff value
__host__ __device__ double Child::getSpeedDiff(const cudaConstants* cConstants) {
    //Calculate the speed diff by calculating the sun-relative speed difference between the spacecraft and the target 
    speedDiff = sqrt(pow(cConstants->vr_fin_target - finalPos.vr, 2) + pow(cConstants->vtheta_fin_target - finalPos.vtheta, 2) + pow(cConstants->vz_fin_target - finalPos.vz, 2)); 

    return speedDiff;
}

// Calculates an orbit posDiff value
// Input: cConstants in accessing properties for the orbit radius of the target
// Output: Assigns and returns this individual's orbitPosDiff value
__host__ __device__ double Child::getOrbitPosDiff(const cudaConstants* cConstants) {
    //Check to see if there is an orbit radius set
    //  Done by seeing if orbital radius is less than 0 since its default value is -1
    if (cConstants->orbitalRadius < 0) {
        //The orbital radius is never set, so set orbitPosDiff to a default bad value
        orbitPosDiff = 1; 
    }
    else {
        //Orbital radius set, calculate the orbital pos diff
        //Calculates the distance between the target and the spacecraft
        double targetToCraftDist = (sqrt(pow(cConstants->r_fin_target, 2) + pow(finalPos.r, 2) + pow(finalPos.z - cConstants->z_fin_target, 2) - (2*finalPos.r*cConstants->r_fin_target*cos(cConstants->theta_fin_target-finalPos.theta))));
        //The distance between the target and the spacecraft should be the orbital radius in order to converge
        orbitPosDiff = abs(targetToCraftDist - cConstants->orbitalRadius);
    }

    return orbitPosDiff;
}

// Calculates an orbit speedDiff value
// Input: cConstants in accessing properties for the orbit speed of the target
// Output: Assigns and returns this individual's orbitSpeedDiff value
__host__ __device__ double Child::getOrbitSpeedDiff(const cudaConstants* cConstants) {
    //Check to see if there is an orbit radius set
    //  Done by seeing if orbital radius is less than 0 since its default value is -1
    if (cConstants->orbitalRadius < 0) {
        //The orbital radius is never set, so set orbitSpeedDiff to a default bad value
        orbitSpeedDiff = 1; 
    }
    else {
        //Orbital radius set, calculate the orbital speed diff
        // The radial velocity of the spacecraft compared with that of the target should be constant to remain in a circular orbit
        //      Thus the pow(cConstants->vr_fin_target - finalPos.vr, 2) should go to 0
        // Additionally, to remain within a circular orbit, the magnitude of the sqrt(v_theta ^2 + v_z^2) must be the orbital velocity
        //      So the (sqrt(v_theta ^2 + v_z^2) - v_orbit)^2 term shoudl also go to 0
        //speedDiff = sqrt(pow(cConstants->vr_fin_target - finalPos.vr, 2) + pow(sqrt(pow(cConstants->vtheta_fin_target - finalPos.vtheta, 2) + pow(cConstants->vz_fin_target - finalPos.vz, 2)) - cConstants->orbitalSpeed, 2));
        orbitSpeedDiff = sqrt(abs((pow(cConstants->vr_fin_target - finalPos.vr, 2) + pow(cConstants->vtheta_fin_target - finalPos.vtheta, 2) + pow(cConstants->vz_fin_target - finalPos.vz, 2)) - pow(cConstants->orbitalSpeed, 2)));
    }

    return orbitSpeedDiff;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

__host__ void Child::getProgress(const cudaConstants* cConstants){

    //Check to see if the individual has an error or not
    //  Previously, if an individual had an error, it would be assigned a progress of 1, which would throw off reporting
    if (errorStatus == VALID || errorStatus == DUPLICATE){
        //The individual is valid, so calcualte the progress normally

        //Holds the progress of the child before it is actually assigned to the child
        //  The relative costs of each objective will be added to it
        //  It will then be divided by the number of objectives and assigned to the child
        double calcProgress = 0; 

        //Iterate through the objectives
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            //Check to see if the goal for this objective is to minimize or maximize the parameter 
            //  Necessary because the progress values are calculated differently depending on the direction 
            //  See the objective header for details on how objective direction is determined
            // if (cConstants->missionObjectives[i].goal < 0) {//Minimization
                
                //See if the child has met the the convergence threshold for this parameter
                if (getParameters(cConstants->missionObjectives[i]) < cConstants->missionObjectives[i].convergenceThreshold) {
                    //Add one to the progress to signify that the parameter has met the goal
                    calcProgress += 1; 
                }
                //The child hasn't met the parameter goal
                else {
                    calcProgress += (getParameters(cConstants->missionObjectives[i])/cConstants->missionObjectives[i].convergenceThreshold); 
                    //Add the progress for this parameter to the goal
                    // if (cConstants->missionObjectives[i].goal < 0) {
                    //     //For minimization, the progress is the parameter divided by the threshold
                    //     calcProgress += (getParameters(cConstants->missionObjectives[i])/cConstants->missionObjectives[i].convergenceThreshold); 
                    // }
                    // else {
                    //     calcProgress += (cConstants->missionObjectives[i].convergenceThreshold/getParameters(cConstants->missionObjectives[i])); 
                    // }
                }
        }

        //The total cost has been calculated
        //It needs to be divided by the number of objectives to find the weighted average progress for each objective
        calcProgress = cConstants->missionObjectives.size()/calcProgress;

        //TODO: the below commented code calculates a progress measure that is easier to understand at a glance, but doesn't work well with our current annealing system. Consider changing progress and annealing in the future

        //Holds the progress of a single objective before it is added to the combined calcProgress variable
        // double objProgress;

        // //Iterate through the objectives
        // for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //     //See if the child has met the the convergence threshold for this parameter
        //     if (getParameters(cConstants->missionObjectives[i]) < cConstants->missionObjectives[i].convergenceThreshold) {
        //         //Set progress to one to signify that the parameter has met the goal
        //         objProgress = 1;
        //     }
        //     //The child hasn't met the parameter goal
        //     else {
        //         //Add the progress for this parameter to the goal
        //         //The progress is calculated by dividing the threshold by the child's parameter value
        //         objProgress = abs((cConstants->missionObjectives[i].convergenceThreshold/getParameters(cConstants->missionObjectives[i]))); 
        //     }

        //     //If the objective is a maximization, store the inverse of the progress for each individual so it is a 0 to 1 scale
        //     if (cConstants->missionObjectives[i].goal > 0){
        //        objProgress = (1/objProgress); 
        //     }
        //     //If the goal is a minimization, the progress for the objective is already on a 0 to 1 scale
            
        //     //Add the objective's progress to the total progress
        //     calcProgress += objProgress;
        // }

        // //The total progress has been calculated
        // //It needs to be divided by the number of objectives to find the weighted average progress for each objective
        // calcProgress /= cConstants->missionObjectives.size();

        //Assign the weighted progress to the child
        progress = calcProgress;
    }
    else {
        //Set the progress to 0 if the individual has an error
        progress = 0; 
    }
}

//NOTE: Long term we need to consider the velocity vector, not just the magnitude
//      this would be where we want to introduce those calculations
//      Three parameters:
//              position diff
//              vector direction (unit vector) diff (this is what we will need to implement)
//              magnitude of velocity diff (this is what we have)
