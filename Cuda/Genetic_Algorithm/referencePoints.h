#ifndef REFERENCEPOINTS_H
#define REFERENCEPOINTS_H

#include "../Config_Constants/constants.h"
#include "adult.h"


///////////////////////////////////////////////////////////////
// ReferencePoints Class & Functions                         //
///////////////////////////////////////////////////////////////

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

    //Constructor which will calculate the values of all of the reference points based on the number of objectives and divisions per objectives
    ReferencePoints(const cudaConstants* cConstants);

    ReferencePoints();

    //Assist function to the constructor which will add a new point to the points vector
    void addPoint(const std::vector<double> values, const cudaConstants* cConstants);
};

//Calculates & returns the normal vector for a plane used to calculate the normalization for any number of objectives greater than or equal to two
std::vector<double> calcNormVector (std::vector<std::vector<double>> matrix, bool initialMatrix);

//Calculates the objective-by-objective relative cost for the adults passed into the function
void calculateRelCost (const cudaConstants *cConstants, ReferencePoints & refPoints, std::vector<Adult> & allAdults);

//Finds the distance between the submitted point and the normalized adult
double findPointDist (const std::vector<double> & point, const Adult & adult);

//Will find the closest reference points to each adult in the newAdults vector
void findAssociatedPoints (const cudaConstants *cConstants, const ReferencePoints & refPoints, std::vector<Adult> & newAdults);

//Function which will handle finding and assigning the rarity of each adult to
void calculateRarity (const cudaConstants *cConstants, const ReferencePoints & refPoints, std::vector<Adult> & allAdults);

#include "referencePoints.cpp"

#endif