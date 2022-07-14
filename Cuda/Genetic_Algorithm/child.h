#ifndef CHILD_H
#define CHILD_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"
#include "../Config_Constants/constants.h"

//TODO: Big step: think about where in the code we rely on numeric data to determine errors (like checking for NAN)
//TODO: Is ERROR_STATUS necessary? We have statuses for it, but they aren't used... we only use STATUS as of 6/28/22

///////////////////////////////////////////////////////////////
// Child Class & Functions                                   //
///////////////////////////////////////////////////////////////

// --- This file includes the definition for the child class and associated functions ---

// Child is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
//once a child is created it will then be copied to an adult with some added parameters
struct Child {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // in AU, difference in position between spacecraft and center of asteroid at end of run
    double speedDiff; // in AU/s, difference in velocity between spacecraft and asteroid at end of run
    double fuelSpent; // in kg, the amount of fuel spent by the individual during its simulation

    double progress; //progress of the individual's posDiff and speedDiff. 0 to 1 scale with 0 being poor and 1 meaning that the individual have completed all objectives
    double avgParentProgress; //The average of the two parents progress

    int birthday; //keeps track of the generation this individual was created in 

    //Both status and error_status are defined in constants.h

    STATUS funcStatus;//flag that keeps track of the status of the child (and eventually adult)
    
    ERROR_STATUS errorStatus; //record of if child is computed correctly, should be set in callRK


//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    // Input: None
    // Outputs: A child with a funcStatus of DEFAULT_CHILD and an errorStatus of NOT_RUN
    // NOTE: Ideally, this should not be used. It serves as a guardrail if for some reason there is a child generated without parameters
    Child();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Called from: generateNewChild()
    // Input:  childParameters - the starting paramters that are generated randomly, from mutation, or from crossover
    //         cConstants - to access the v_escape value for calculations of the startParams
    //         genCreated - the generation this child is created, used for birthday/age data
    //         calcAvgParentProgress - the creating parent's average progress. Will be set to 0 for randomly generated children
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    //         Child is ready to be passed into callRK
    Child(rkParameters<double> & childParameters, const cudaConstants* cConstants, int genCreated, double calcAvgParentProgress);

    // Copy constructor
    // Sets this child's parameters, elements, posDiff, and speedDiff to another child's values for these quantities
    // Input: Another child
    // Output: this child is a copy of the other child that was passed in
    //TODO: See if we need this defined (the default might work?)
    Child(const Child& other);

    #ifdef UNITTEST //set up in unit_testing_main.cpp
    // Child constructor only with status - ONLY USED FOR UNIT TESTING!!!
    // Input: errorStatus - this is the only thing set in Child by this constructor
    // Output: A child with only its errorStatus set
    // ONLY USED FOR UNIT TESTING!!! DO NOT USE ELSEWHERE
    Child(ERROR_STATUS status): errorStatus(status){}

    // Child constructor for testing_genetics so we don't have to callRK on a child
    // Input: errorStatus, posDiff, and speedDiff
    // Output: a child that can be converted into an adult and be given a rank and distance without going through callRK
    Child(double posD, double speedD): errorStatus(VALID), posDiff(posD), speedDiff(speedD), funcStatus(FUNCTIONAL_CHILD){}

    // Child constructor for testing_genetics that doesn't take cConstants as a parameter because we're not using cConstants in testing_genetics
    // If I need to have elements for some reason, this is easy to modify to shadow the other Child constructor with rkParameters, I just need the actual v_escape value
    Child(rkParameters<double> & childParameters): startParams(childParameters){}

    // Child constructor for unit tests that makes it possible to only use tripTime, posDiff, and speedDiff as input parameters
    //alpha is 3.12, beta is 1.5, and zeta is 0 -> these were all arbitrarily chosen values
    Child(double tripTime, double posD, double speedD): startParams(tripTime, 3.12, 1.5,0, coefficients<double>()), errorStatus(VALID), posDiff(posD), speedDiff(speedD), funcStatus(FUNCTIONAL_CHILD){}
   
    // Child constructor for testing_genetics that will copy all the important elements from an Adult back into a Child
    //This allows us to use the unit testing verification functions to ensure newAdults were correctly generated, even though those functions are meant to take in Children
    Child(rkParameters<double> & childParameters, double posD, double speedD, double prog, double avgPrntProg, int bday): startParams(childParameters), errorStatus(VALID), posDiff(posD), speedDiff(speedD), funcStatus(FUNCTIONAL_CHILD), progress(prog), avgParentProgress(avgPrntProg), birthday(bday){}

    #endif //unit testing end if

    //Child(rkParameters<double> & childParams, elements<double> posFinal,  double posD,  double speedD, STATUS s, int errStat): startParams(childParams), finalPos(posFinal), posDiff(posD), speedDiff(speedD), funcStatus(s), errorStatus(errorStatus){}

    rkParameters<double> getRKParams() {return startParams;}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    //Getter for a parameter dependent on the objective that is passed in
    //It will return the parameter needed for comparisons and the eventual sorting of adults
    //  Parameters include diffs, fuel spent, and triptime
    // Inputs: requestObjective - the objective which will be used to grab the right parameter from the child
    // Output: the parameter that is relevent to the objective
    //              Could be posDiff, fuelSpent, tripTime, etc
    __host__ double getParameters (const objective & requestObjective) const;

    // Calculates a posDiff value
    // Input: cConstants in accessing properties for the final position of the target (such as r_fin_ast, theta_fin_ast, and z_fin_ast)
    // Output: Assigns and returns this individual's posDiff value
    __host__ __device__ double getPosDiff(const cudaConstants* cConstants);

    // Calculates a speedDiff value
    // Input: cConstants in accessing properties for the final velocity of the target (such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast)
    // Output: Assigns and returns this individual's speedDiff value
    __host__ __device__ double getSpeedDiff(const cudaConstants* cConstants);

    // Calculates the progress depending on the mission type
    // Input: cConstants for the tolerances, and the child's posDiff and speedDiff are used
    //        NOTE: This function assumes that the pos/speed diffs have already been calculated
    // Output: the child has been assigned a progress value betweeen 0 and 1
    __host__ void getProgress(const cudaConstants* cConstants);
};

#include "child.cpp"

#endif
