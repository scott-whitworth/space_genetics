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

    bool isClone; //bool that tells if the individual is a clone of the best cost individual or not.
    double difference; //represents the percent similar to the best cost individual

    bool isParent; //Indicates if this individual is already chosen as a survivor

    int dominatedCount;
    int rank;
    std::vector<int> dominated;
    double distance;

    // Default constructor
    Individual();

    // Set the initial position of the spacecraft according to the newly generated parameters
    // Input: cConstants - to access c3energy value used in getCost()
    //        newInd - struct returned by generateNewIndividual()
    // Output: this individual's startParams.y0 is set to the initial position and velocity of the spacecraft
    Individual(rkParameters<double> & newInd, const cudaConstants* cConstants);

    // Calculates a posDiff value
    // Input: cConstants in accessing properties such as r_fin_ast, theta_fin_ast, and z_fin_ast
    // Output: Assigns and returns this individual's posDiff value
    __host__ __device__ double getPosDiff(const cudaConstants* cConstants);

    // Calculates a speedDiff value
    // Input: cConstants in accessing properties such as vr_fin_ast, vtheta_fin_ast, and vz_fin_ast
    // Output: Assigns and returns this individual's speedDiff value
    __host__ __device__ double getSpeedDiff(const cudaConstants* cConstants);

    // Calculates a cost value to quantitatively evaluate this Individual
    // Input: cConstants in accessing properties such as pos_threshold, c3energy, and v_impact
    // Output: Assigns and returns this individual's cost value
    //__host__ __device__ double getCost(const cudaConstants* cConstants);
    __host__ __device__ double getCost_Hard(const cudaConstants* cConstants);
    __host__ __device__ double getCost_Soft(const cudaConstants* cConstants);

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

//bool BetterSpeedDiff(Individual& personA, Individual& personB);
bool HigherSpeedDiff(Individual& personA, Individual& personB);
bool LowerSpeedDiff(Individual& personA, Individual& personB);

//Utility to calculate the position difference between a position velocity set and the asteroid final position
// Input: currentState - set of position and velocity coordinates
//                       Often the final position/velocity of the spacecraft at the end of a RK run
//        cConstants   - Configuration information, primarily using the Asteroid final position and velocity
// Return: Root of the sum of squares
//         posDiff = sqrt(  (ast_r - craft_r) ^ 2 + (ast_r * ast_theta - craft_r * craft_theta % 2Pi)^2 +  (ast_z - craft_z)^2  )
//__host__ __device__ double calcPosDiff(const elements<double>& currentState, const cudaConstants* cConstants);

double checkDifference(Individual& personA, Individual& personB);

bool same(Individual& personA, Individual& personB, double percent);

bool notAClone(Individual& personA, Individual& personB);

bool dominates(Individual& personA, Individual& personB);

bool rankSort(Individual& personA, Individual& personB);

#include "individuals.cpp"

#endif
