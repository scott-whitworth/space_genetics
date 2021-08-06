#ifndef INDIVIDUALS_H
#define INDIVIDUALS_H

#include "../Runge_Kutta/rkParameters.h"
#include "../Earth_calculations/earthInfo.h"

// Individual is a structure member of the genetic algorithm's population and has set of parameters and the resulting position and velocity
struct Individual {
    rkParameters<double> startParams; // input parameters for the run- the unique identifiers of this individual

    elements<double> finalPos; // final position of the spacecraft at end of run

    double posDiff; // difference in position between spacecraft and center of asteroid at end of run
    double speedDiff; // difference in velocity between spacecraft and asteroid at end of run
    double cost;    // cost value of the individual, something that the genetic algorithm is attempting to minimize
    bool isParent; //Indicates if this individual is already chosen as a survivor

    int dominatedCount; //keeps track of how many individuals this individual has been dominated by
    int rank; //how good this individual is compared to other individuals in it's generation. The lower the rank, the better.
    std::vector<int> dominated; //The pool indexes of individuals that have been dominated by this individual 
    double distance; //how different this individual is from the individuals most similar to it.

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Default constructor
    Individual();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Input: cConstants - to access c3energy value used in getCost()
    //        newInd - struct returned by generateNewIndividual()
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    Individual(rkParameters<double> & newInd, const cudaConstants* cConstants);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Calculates a posDiff value
    // Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
    // Output: Assigns and returns this individual's posDiff value
    __host__ __device__ double getPosDiff(const cudaConstants* cConstants);

    // Calculates a speedDiff value
    // Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
    // Output: Assigns and returns this individual's speedDiff value
    __host__ __device__ double getSpeedDiff(const cudaConstants* cConstants);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Calculates a cost value to quantitatively evaluate this Individual
    // Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
    // Output: Assigns and returns this individual's cost value
    __host__ __device__ double getCost_Hard(const cudaConstants* cConstants);
    __host__ __device__ double getCost_Soft(const cudaConstants* cConstants);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Comparison operators overloaded to compare individuals by their cost values (determined in getCost())
    bool operator>(Individual &other);
    bool operator<(Individual &other);
    bool operator==(Individual &other);
};

// Compare two individuals by their positional difference values
// input: two individuals
// output: returns true if personB has a higher positional difference than personA
bool LowerPosDiff(Individual& personA, Individual& personB);

// Compare two individuals by their velocity difference values
// input: two individuals
// output: returns true if personB has a higher velocity difference than personA
bool HigherSpeedDiff(Individual& personA, Individual& personB);

// Compare two individuals by their velocity difference values
// input: two individuals
// output: returns true if personB has a lower velocity difference than personA
bool LowerSpeedDiff(Individual& personA, Individual& personB);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals to see if the first individual dominates the second individual
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
bool dominates(Individual& personA, Individual& personB, const cudaConstants* cConstants);

//Compare two individuals by their rank
//WARNING: Using this function to sort individuals will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the individuals are rank 1, sorting them using this method will do nothing. 
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(Individual& personA, Individual& personB);

//Compare two individuals by their rank and distance
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a greater distance than person B, return true
bool rankDistanceSort(Individual& personA, Individual& personB);

#include "individuals.cpp"

#endif
