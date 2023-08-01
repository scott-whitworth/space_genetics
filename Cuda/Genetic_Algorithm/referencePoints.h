#ifndef REFERENCEPOINTS_H
#define REFERENCEPOINTS_H

#include "../Config_Constants/constants.h"
#include "adult.h"
#include "../Genetic_Algorithm/sort.h"


///////////////////////////////////////////////////////////////
// ReferencePoints Class & Functions                         //
///////////////////////////////////////////////////////////////

// Calculations based off of NGSAIII algorithm

// --- This file includes the class & functions which handle reference points, including the calculations of the points and assigning rarity based on the points ---

//ReferencePoints is a class which handles the creation of reference points and calculations relating to finding the rarity for all individuals
class ReferencePoints {

    public:

    //2D vector which will store all of a points and their individual components
    //  The first dimension will store each point
    //  The second dimension will store the values for each objective
    std::vector<std::vector<double>> points;

    //Vectors that will assist with normalization
    //  To do proper normalization, both the overall worst and best value for each objective need to be stored
    //Store the adults with the worst values
    std::vector<Adult> objWorst;
    //Store the adults with the best values
    std::vector<Adult> objBest;

    //!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Constructor which will calculate the values of all of the reference points based on the number of objectives and divisions per objectives
    // Input: cConstants - used to access the number of divisions and the number of objectives 
    // Output: Fills the points vector with all calculated reference points   
    ReferencePoints(const cudaConstants* cConstants);

    //Default constructor that should not be used
    ReferencePoints();

    //!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Functions involved with calculating the reference points

    //Assist function to the constructor which will add a new point to the points vector
    // Input:   values - holds the coordinate values for the reference point that will be created 
    //          cConstants - used to access the number of divisions and the number of objectives 
    // Output: adds a new point to the points vector based on the passed in values  
    void addPoint(const std::vector<double> values, const cudaConstants* cConstants);
};

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Functions involved with normalizing the adult's objective values

//Calculates & returns the normal vector for a plane used to calculate the normalization for any number of objectives greater than or equal to two
// Input:   matrix - 2D vector passed in that is used for the calculation 
//          initialMatrix - determines if the function call is the first call (return the normal vector coordinates) or not (return the determinant of the vector) 
// Output: Overall it will return the normal vector for the passed in matrix, but it will also return the determinant for the matrix to itself (recursive function)
std::vector<double> calcNormVector (std::vector<std::vector<double>> matrix, bool initialMatrix);

//Calculates the objective-by-objective relative cost for the adults passed into the function
// Input:   cConstants - used to access the number of objectives and thresholds for the objectives
//          refPoints - the algorithm's reference points
//          allAdults - the vector of adults which will have relative costs calculated for them
// Output: Each adult will have relative costs for each objective calculated and assigned to them   
void calculateRelCost (const cudaConstants *cConstants, ReferencePoints & refPoints, std::vector<Adult> & allAdults);

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Functions involved with calculating the rarity for each adult

//Finds the distance between the submitted point and the normalized adult
// Input:   point - the reference point to be compared against the adult
//          adult - the adult of interest and its relative costs
// Output: returns the distance between the relative costs of the adult and the line between the origin and the reference point  
double findPointDist (const std::vector<double> & point, const Adult & adult);

//Will find the closest reference points to each adult in the newAdults vector
// Input:   refPoints - the list of reference points 
//          newAdults - the list of adults which will have the closest reference point found
// Output: each adult will be assigned with the index of the reference point within refPoints which is closest to it  
void findAssociatedPoints (const ReferencePoints & refPoints, std::vector<Adult> & newAdults);

//Function which will handle finding and assigning the rarity of each adult to
// Input:   refPoints - the list of reference points 
//          newAdults - which list of adults which will be compared to find relative rarity
// Output: each adult will be assigned with a rarity score depending on which   
void calculateRarity (const ReferencePoints & refPoints, std::vector<Adult> & allAdults);

#include "referencePoints.cpp"

#endif