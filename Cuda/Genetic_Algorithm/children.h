#ifndef CHILDREN_H
#define CHILDREN_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"

//enumeration to make error status easier to keep track of, as opposed to using hard-coded numbers 
enum STATUS {  //TODO: Possibly move to constants.h
    VALID = 0, //TODO: Document these
    SUN_ERROR = 1,
    OTHER_ERROR = 2,
};

//TODO: Big step: think about where in the code we rely on numeric data to determine errors (like checking for NAN)
//      Change to interface with STATUS instead

//Used to fill bad values in the pos and speed diffs of individuals with errors
enum ERROR_VALUES { //TODO: these should be #define in config.h
    BAD_POSDIFF = 10,
    BAD_SOFT_SPEEDDIFF = 10, //TODO: I might change to impact / rendev.
    BAD_HARD_SPEEDDIFF = 0,
};

//TODO: Clarify this header
// Individual is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
struct Child {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // in AU, difference in position between spacecraft and center of asteroid at end of run
    double speedDiff; // in AU/s, difference in velocity between spacecraft and asteroid at end of run

    //flag that keeps track of the status of the child (and eventually adult)
    STATUS status;

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    //TODO: what does the default do?
    Child();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Called from: generateNewChild()
    // Input: cConstants - to access c3energy value used in getCost()
    //TODO: Confirm that we need cConstants in this call (if not, remove)
    //        newChild - struct returned by generateNewChild()
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    //         Child is ready to be passed into callRK
    //TODO: newChild probably needs a new name (like childParameters)
    Child(rkParameters<double> & newChild, const cudaConstants* cConstants);

    // Copy constructor
    // Sets this child's parameters, elements, posDiff, and speedDiff to another child's values for these quantities
    // Input: Another child
    // Output: this child is a copy of the other child that was passed in
    //TODO: See if we need this defined (the default might work?)
    Child(const Child& other);

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

#include "children.cpp"

#endif
