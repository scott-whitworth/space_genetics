#ifndef CHILD_H
#define CHILD_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"

//enumeration to make error status easier to keep track of, as opposed to using hard-coded numbers 
//could be removed in the future, included for now to check if the default constructor is being used
enum STATUS {
    DEFAULT_CHILD = 0, //child that is created through the default constructor, not ready to be made an adult
    FUNCTIONAL_CHILD = 1, //child that is created in main constructor and given startParams
    FUNCTIONAL_ADULT = 2, //once the child has gone through callRK and become an adult, we can change its status to FUNCTIONAL_ADULT
};

//TODO: Big step: think about where in the code we rely on numeric data to determine errors (like checking for NAN)
//      Change to interface with STATUS instead

//TODO: Clarify this header
// Individual is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
struct Child {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // in AU, difference in position between spacecraft and center of asteroid at end of run
    double speedDiff; // in AU/s, difference in velocity between spacecraft and asteroid at end of run

    STATUS funcStatus;//flag that keeps track of the status of the child (and eventually adult)
    
    int errorStatus;//flag that keeps track of what kind of nan it is and if it is a nan

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    //TODO: what does the default do?
    Child();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Called from: generateNewChild()
    // Input: cConstants - to access the v_escape value for calculations of the startParams
    //        newChild - struct returned by generateNewChild()
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    //         Child is ready to be passed into callRK
    //TODO: newChild probably needs a new name (like childParameters)
    Child(rkParameters<double> & childParameters, const cudaConstants* cConstants);

    // Copy constructor
    // Sets this child's parameters, elements, posDiff, and speedDiff to another child's values for these quantities
    // Input: Another child
    // Output: this child is a copy of the other child that was passed in
    //TODO: See if we need this defined (the default might work?)
    Child(const Child& other);

    // Child constructor only with status - ONLY USED FOR UNIT TESTING!!!
    // Input: status - this is the only thing set in Child by this constructor
    // Output: A child with only its errorStatus set
    // ONLY USED FOR UNIT TESTING!!! DO NOT USE ELSEWHERE
    Child(int status): errorStatus(status){}

    //Child(rkParameters<double> & childParams, elements<double> posFinal,  double posD,  double speedD, STATUS s, int errStat): startParams(childParams), finalPos(posFinal), posDiff(posD), speedDiff(speedD), funcStatus(s), errorStatus(errorStatus){}



//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    // Calculates a posDiff value
    // Input: cConstants in accessing properties for the final position of the target (such as r_fin_ast, theta_fin_ast, and z_fin_ast)
    // Output: Assigns and returns this individual's posDiff value
    __host__ __device__ double getPosDiff(const cudaConstants* cConstants);

    // Calculates a speedDiff value
    // Input: cConstants in accessing properties for the final velocity of the target (such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast)
    // Output: Assigns and returns this individual's speedDiff value
    __host__ __device__ double getSpeedDiff(const cudaConstants* cConstants);

};

#include "child.cpp"

#endif
