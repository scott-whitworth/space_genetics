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
    //Worst objective values
    std::vector<double> objWorst;
    //Best objective values
    std::vector<double> objBest;

    //Constructor which will calculate the values of all of the reference points based on the number of objectives and divisions per objectives
    ReferencePoints(const cudaConstants* cConstants);

    ReferencePoints();

    //Assist function to the constructor which will add a new point to the points vector
    void addPoint(const std::vector<double> values, const cudaConstants* cConstants);
};

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