#include <math.h>
#define _USE_MATH_DEFINES //For INT_MAX

//Constructor which will calculate the values of all of the reference points based on the number of objectives and divisions per objectives
ReferencePoints::ReferencePoints(const cudaConstants* cConstants) {
    //Create the vector that will be used to calculate the new values
    std::vector<double> newValues(cConstants->missionObjectives.size() - 1, 0);

    //Tracking variables
    int valTot = 0;      //Stores the total value inside the newValues vector so the total value doesn't go too high
    int indexTrack = newValues.size() - 1; //Index tracking cursor to move through the newValues vector when modifying its values

    //Create new points while the first value in the new point is less than the number of divisions
    while (newValues[0] < cConstants->divisions) {

        //Create a new point based on the current values
        addPoint(newValues, cConstants);

        //Add one to the last index in the new values array to start creating the next point
        newValues[newValues.size()-1] += 1;
        valTot += 1;

        //While loop to make sure the total value of the vector doesn't exceed the maximum value (based on cConstant's divisions)
        while (valTot > cConstants->divisions) {
            //Reset the value of the tracked index
            valTot -= newValues[indexTrack];
            newValues[indexTrack] = 0;

            //Go to the previous index and add one to its value
            indexTrack--;
            newValues[indexTrack] += 1;
            valTot += 1;
        }

        //Reset the index cursor
        indexTrack = newValues.size() - 1;
    }

    //Add the last point (where the first value is 1 and the rest are 0)
    addPoint(newValues, cConstants);
}

void ReferencePoints::addPoint(const std::vector<double> values, const cudaConstants* cConstants){
    //Calculate what the value of the last point needs to be
    //  Do this by calculating the value of 1 - the combined value of the rest of the values vector
    double valTot = 0;

    //New point vector
    std::vector<double> newPoint;

    //Calculate the total value of values and push the values to the new point
    //All values in the new point need to be divided by the number of divisions to normalize the values to a 0 to 1 scale
    for (int i = 0; i < values.size(); i++){
        valTot += values[i];
        newPoint.push_back((values[i])/cConstants->divisions);
    }

    //Add the last value the new point
    newPoint.push_back((cConstants->divisions-valTot)/cConstants->divisions);

    //Push the new point to the points vector
    points.push_back(newPoint);

    return;
}

//Calculates & returns the normal vector for a plane used to calculate the normalization for any number of objectives greater than or equal to two
std::vector<double> calcNormVector (std::vector<std::vector<double>> matrix, bool initialMatrix) {

    //vector that will either hold the normal vector, or just the determinant value
    //  if initialMatrix is true, this will hold the normal vector
    //  if initialMatrix is fale, this will hold the determinant value
    std::vector<double> norm;

    //If the matrix size of 2, the recursion ends
    //  The determinant of the matrix can be directly calculated
    //  Or the mission has two objectives, so finding the normalized vector is trivial
    if (matrix.size() == 2) { 

        //Mission has two objectives, add the components of the normal vector
        if (initialMatrix) { 
            norm.push_back(matrix[1][1]); 
            norm.push_back(-matrix[1][0]);
        } 
        //the matrix is 2x2, but this is a recursive call, calculate the determinant
        else { 
            norm.push_back(((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]))); 
        } 

        //Return the normal or the determinant
        return norm; 
    } 

    //Need to make a recursive call to find the full normal
    else { 

        //For each value across the top row of the matrix 
        for (int i = 0; i < matrix.size(); i++) { 

            //Reset the determinant to be a size n-1 square matrix
            std::vector<std::vector<double>> det (matrix[0].size()-1, std::vector<double>(matrix.size()-1, 0)); 

            //Initialize column and row trackers for the determinant
            int col = 0, row = 0; 

            //Add the value at column j (not including the top row)... 
            for (int j = 0; j < matrix.size(); j++) { 

                //only add if the column is not the same one as the top row value 
                if (j != i) { 

                    //...and row k... 
                    for (int k = 1; k < matrix.size(); k++) { 

                        //...into the determinant matrix
                        det[col][row] = matrix[k][j]; 

                        //Add the next value to the next column
                        col++; 

                    }

                    //add the new values to the next row
                    row ++; 
                    //Reset the column counter to start at the beginning of the row 
                    col = 0; 

                } 

            } 

            //Calculate the determinant 
            std::vector <double> calcDet = calcNormVector(det, false); 

             
            //If this isn't the first time calling the function, norm will store the determinant, not the normalized vector
            if (!initialMatrix) { 
                //Make sure there is a space to add the determinant
                norm.resize(calcDet.size()); 
            }

            //Determine whether to add the positive or negative value of the determinant 
            if (i % 2 == 0) { 
                //Add the determinant

                //If this is the initial call, the whole vector needs to be stored instead of calculating a determinant, so push the calculated determinant to the back of the norm vector
                if (initialMatrix) { 
                    norm.push_back(calcDet[0]); 
                } 

                //Not the first call, calculate this part of the determinant
                else { 
                    norm[0] += (matrix[0][i]*calcDet[0]); 
                } 
            } 

            else { 
                //Subtract the determinant

                //If this is the initial call, the whole vector needs to be stored instead of calculating a determinant, so push the calculated determinant to the back of the norm vector
                if (initialMatrix) { 
                   norm.push_back(-(calcDet[0]));  
                } 

                //Not the first call, calculate this part of the determinant
                else { 
                    norm[0] += -(matrix[0][i]*calcDet[0]); 
                } 
            } 
        } 
    } 

    //Return the normal vector or the determinant
    return norm;
}

//Calculates the objective-by-objective relative cost for the adults passed into the function
void calculateRelCost (const cudaConstants *cConstants, ReferencePoints & refPoints, std::vector<Adult> & allAdults) {

    //See if the adult vectors need to be filled in
    //  If the size of the best/worst value vectors is 0, it means this is the first generation and they need to be set to an arbitrary adult
    // if (refPoints.objBest.size() == 0) {
        //Make sure the vectors are cleared
        refPoints.objWorst.clear();
        refPoints.objBest.clear();

        //Push back a random adult for each objective
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            refPoints.objWorst.push_back(allAdults[0]);
            refPoints.objBest.push_back(allAdults[0]);
        }
    // }

    //Find the normalizations for each objective
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //Go through the rest of the adults to see if there are new best/worst values
        for (int j = 0; j < allAdults.size(); j++) {

            // std::cout << "\nTEST: comparing value of " << allAdults[j].getParameters(cConstants->missionObjectives[i]) << " for objective "  << i;
            // std::cout << " against worst " << refPoints.objWorstValue[i];
            // std::cout << " and best " << refPoints.objBestValue[i] << ".\n";

            //Check to see if the adult has the new worst value for this objective
            //Because maximizations are stored as negatives, worse values are going to be larger values, regardless of objective
            if (allAdults[j].getParameters(cConstants->missionObjectives[i]) > refPoints.objWorst[i].getParameters(cConstants->missionObjectives[i])) {

                // std::cout << "\nTEST: Found new worst value " << allAdults[j].getParameters(cConstants->missionObjectives[i]) << " of obj " << i << " at adult " << j << ".\n";

                //Found a new worse value, store the adult
                refPoints.objWorst[i] = allAdults[j];
            }
            //Check if the adult has the new best value for this objective
            else if (allAdults[j].getParameters(cConstants->missionObjectives[i]) < refPoints.objBest[i].getParameters(cConstants->missionObjectives[i])) {

                // std::cout << "\nTEST: Found new best value " << allAdults[j].getParameters(cConstants->missionObjectives[i]) << " of obj " << i << " at adult " << j << ".\n";

                //Found a new best value, store the adult
                refPoints.objBest[i] = allAdults[j];
            }          
        }
    }


    //Holds the points for the objective intercepts
    std::vector<double> intercepts;

    //Find intercepts if there are three objectives
    if (cConstants->missionObjectives.size() >= 2){
        // std::cout << "\nUsing plane intercepts:\n";

        // //First dimension is the individual (1,2,3)
        // //Second dimension is the coordinate (x,y,z)
        // std::vector<std::vector<double>> worstPoints;

        // for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //     //Push back a new set of worst points for this objective
        //     worstPoints.push_back(std::vector<double>(0,0));

        //     //For the objective's worst individual, store its objective-specific values
        //     for (int j = 0; j < cConstants->missionObjectives.size(); j++){
        //         if (cConstants->missionObjectives[j].goal < 0) {
        //             worstPoints[i].push_back(refPoints.objWorst[i].getParameters(cConstants->missionObjectives[j]));
        //         } 
        //         else {
        //             worstPoints[i].push_back(-refPoints.objWorst[i].getParameters(cConstants->missionObjectives[j]));
        //         }

        //     std::cout << "\tPoint " << i << " " << j << ": " << worstPoints[i][j] << "\n";
        //     }
        // }

        // // Coefficients for the unit vector for each objective
        // double c1 = (((worstPoints[1][1]-worstPoints[0][1]) * (worstPoints[2][2]-worstPoints[0][2])) - ((worstPoints[2][1]-worstPoints[0][1]) * (worstPoints[1][2]-worstPoints[0][2]))); //scalar in i direction (y2-y1)(z3-z1) - (y3-y1)(z2-z1)
        // double c2 = (((worstPoints[2][0]-worstPoints[0][0]) * (worstPoints[1][2]-worstPoints[0][2])) - ((worstPoints[1][0]-worstPoints[0][0]) * (worstPoints[2][2]-worstPoints[0][2]))); //scalar in j direction (x3-x1)(z2-z1) - (x2-x1)(z3-z1)
        // double c3 = (((worstPoints[1][0]-worstPoints[0][0]) * (worstPoints[2][1]-worstPoints[0][1])) - ((worstPoints[2][0]-worstPoints[0][0]) * (worstPoints[1][1]-worstPoints[0][1]))); //scalar in k direction (x2-x1)(y3-y1) - (x3-x1)(y2-y1)

        // //Check to make sure that the scalars are not zero
        // // if (abs(c1) < cConstants->doublePrecThresh) {
        // //     c1 = cConstants->doublePrecThresh;
        // // }
        // // if (abs(c2) < cConstants->doublePrecThresh) {
        // //     c2 = cConstants->doublePrecThresh;
        // // }
        // // if (abs(c3) < cConstants->doublePrecThresh) {
        // //     c3 = cConstants->doublePrecThresh;
        // // }

        // std::cout << "\n\tScalar 1: " << c1 << "\n";
        // std::cout << "\tScalar 2: " << c2 << "\n";
        // std::cout << "\tScalar 3: " << c3 << "\n";

        // //Calculate the objective intercepts
        // intercepts[0] = (((c1*worstPoints[0][0]) + (c2*worstPoints[0][1]) + (c3*worstPoints[0][2]))/c1); //x intercept, (c1x1 + c2y1 + c3z1)/c1
        // std::cout << "\n\tIntercept 1: " << intercepts[0] << "\n";

        // intercepts[1] = (((c1*worstPoints[0][0]) + (c2*worstPoints[0][1]) + (c3*worstPoints[0][2]))/c2); //y intercept, (c1x1 + c2y1 + c3z1)/c2
        // std::cout << "\tIntercept 2: " << intercepts[1] << "\n";

        // intercepts[2] = (((c1*worstPoints[0][0]) + (c2*worstPoints[0][1]) + (c3*worstPoints[0][2]))/c3); //z intercept, (c1x1 + c2y1 + c3z1)/c3
        // std::cout << "\tIntercept 3: " << intercepts[2] << "\n\n";




        //Get the on-plane vectors
        //Create a matrix of n, n-sized vectors between the worst points
        std::vector<std::vector<double>> matrix;

        std::cout << "\n\nTEST - created matrix:";

        //Calculate the on-plane vectors
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            //Add a new empty vector to the matrix
            matrix.push_back(std::vector<double>(0,0));

            std::cout << "\n";

            for (int j = 0; j < cConstants->missionObjectives.size(); j++) {
                //The vectors are all relative to the first objective's worst adult's points
                //  so each component will be j adult's point minus the 0th adult's point
                //  the top row of the matrix will be all zeroes, but it won't matter for the calculation of the normal vector
                matrix[i].push_back(refPoints.objWorst[i].getParameters(cConstants->missionObjectives[j]) - refPoints.objWorst[0].getParameters(cConstants->missionObjectives[j]));

                std::cout << "  " << matrix[i][j];
            }
        }

        //The matrix is created, calculate the plane's normal vector
        std::vector<double> normal = calcNormVector(matrix, true); 

        std::cout << "\n\nNormal vector: ";

        //Calculate the intercepts
        // The formula for intercept x (where the 0th adult is a, the normal objective is n, and the objectives are 1,2,3...) is (a1n1 + a2n2 + ...)/nx

        //Value to store the numerator for the intercepts
        double num = 0;

        //The numerator is consistent, use a for loop to calculate it
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            std::cout << "\n  " << normal[i];

            num += (refPoints.objWorst[0].getParameters(cConstants->missionObjectives[i]) * normal[i]);
        }

        std::cout << "\n\nIntercepts: ";

        //Calculate the intercepts
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            intercepts.push_back(num / normal[i]);

            std::cout << "\n  " << intercepts[i];
        }

        std::cout << "\n\n";
    }

    //The normalization values have been found, go through the adults and calculate the objective costs
    for (int i = 0; i < allAdults.size(); i++) {

        //Calculate the cost for each objective
        for (int j = 0; j < cConstants->missionObjectives.size(); j++){
            //Establish numerator and denominator for the normalization calculation
            //The numerator is the individual's objective value minus the best found objective value
            double num = allAdults[i].getParameters(cConstants->missionObjectives[j]) - refPoints.objBest[j].getParameters(cConstants->missionObjectives[j]);
            //The denominator is dependent on the configuration
            double denom;

            //If the normalization is based on objective intercepts (happens if there are two or more objectives), denominator is objective intercept - best found objective value
            if (cConstants->missionObjectives.size() >= 3){
                denom = intercepts[j] - refPoints.objBest[j].getParameters(cConstants->missionObjectives[j]);
            }
            //If the normalization is not based on objective intercepts, the denominator is the worst objective value - the best objective value
            else {
                denom = refPoints.objWorst[j].getParameters(cConstants->missionObjectives[j]) - refPoints.objBest[j].getParameters(cConstants->missionObjectives[j]);
            }

            //Make sure there are no divide by 0 errors
            if (abs(denom) < cConstants->doublePrecThresh) {
                //minimization vs maximization handling
                if (cConstants->missionObjectives[j].goal < 0) {
                    denom = cConstants->doublePrecThresh;
                }
                else {
                    denom = -cConstants->doublePrecThresh;
                }
            }

            //Calculate the actual normailized objective
            allAdults[i].normalizedObj[j] = num/denom;
        }
    }

    return;
}

//Finds the distance between the submitted point and the normalized adult
double findPointDist (const std::vector<double> & point, const Adult & adult) {
    //The distance between a point and a line where A = point on line (origin), B = second point on line (reference point), P = point of interest (adult)
    //  pa = P - A; -> pa = P (adult point)
    //  ba = B - A; -> ba = B (reference point)
    //  t = dot(pa, ba)/dot(ba, ba);
    //  your_dist = norm(pa - t * ba,2)

    //holds the dot product of pa & ba
    double numDot = 0;
    //holds the dot product of ba with itself
    double denDot = 0;

    //Calculate the dot products
    for (int i = 0; i < point.size(); i++) {
        numDot += (point[i] * adult.normalizedObj[i]);
        denDot += (point[i] * point[i]);
    }

    //Calculate t with the dot products
    double t = numDot/denDot;

    //Calculate the distance
    double distance = 0;

    //The norm is the square of each of the components squared
    //Calculate the sum of each compenent squared]
    for (int i = 0; i < point.size(); i++) {
        distance += pow(adult.normalizedObj[i] - (t * point[i]), 2);
    }

    //The norm (and thus distance) is the square root of the sum
    distance = sqrt(distance);

    //return the distance
    return distance;
}

//Will find the closest reference points to each adult in the newAdults vector
void findAssociatedPoints (const cudaConstants *cConstants, const ReferencePoints & refPoints, std::vector<Adult> & newAdults) {
    //Trackers for finding the closest reference point
    //In the normalized coordinates, the farthest corner (1,1,1) is at a distance of sqrt(3)=1.732... from the origin.
    double minDistVal = 2;      //Tracks the distance of the reference point currently the closest from the adult

    //Value that holds the distance from the adult to the reference point being checked
    double curDist = 2;

    //Find the closest point for all newAdults
    for (int i = 0; i < newAdults.size(); i++) {

        //For each adult check each reference point
        for (int j = 0; j < refPoints.points.size(); j++) {
            
            //Calculate the distance to the reference point
            //First add the squares of the differences for each objective/component
            // for (int k = 0; k < refPoints.points[0].size(); k++) {
                
            //     curDist += pow((refPoints.points[j][k] - newAdults[i].normalizedObj[k]), 2);
            // }

            //Finally take the square root to find the distance
            // curDist = sqrt(curDist);
            curDist = findPointDist(refPoints.points[j], newAdults[i]);

            //See if it is better than the current minimum distance
            if (curDist <  minDistVal) {
                //Store the new minDistVal
                minDistVal = curDist;

                //Store the new closest point
                newAdults[i].associatedRefPoint = j;
            }
        }
    }

   return; 
}

//Function which will handle finding and assigning the rarity of each adult
void calculateRarity (const cudaConstants *cConstants, const ReferencePoints & refPoints, std::vector<Adult> & allAdults) {
    //2D vector which will keep track of which adults are closest to which reference point
    //  The first dimension corresponds to the indexes of the reference points
    //  The second dimension is a list of the indexes of all adults who have the reference point as closest to them
    std::vector<std::vector<int>> refPointAssociations (refPoints.points.size(), std::vector<int>(0));

    //Vector which will keep track of how many adults are associated with each reference point
    //  The indexes of this vector correspond to the indexes of the reference points
    std::vector<int> refPointAssocNum (refPoints.points.size(), 0);

    //Go through all adults to generated the association lists
    for (int i = 0; i < allAdults.size(); i++) {
        //Add the adult to the list of their closest reference point
        refPointAssociations[allAdults[i].associatedRefPoint].push_back(i);

        //Add one to the point's number of associated adults
        refPointAssocNum[allAdults[i].associatedRefPoint]++;
    }

    //Rarity counter will both be used to assign rarity scores and determine when all rarity scores have been assigned
    int rarityCount = 0;

    //Trackers for assigning rarity
    int minAssocNum = allAdults.size() + 1; //Tracks the minimum number of found associations
    int minAssocIndex = 0;                  //Tracks the index where the minimum number of associations was found

    //Loop while there are adults with unassigned rarity values
    while (rarityCount < allAdults.size()) {
        //Go through all points and find which one has the fewest associations
        for (int i = 0; i < refPointAssocNum.size(); i++) {

            //Check to see if the current reference point has a new minimum amount of associations
            if (refPointAssocNum[i] < minAssocNum) {
                //Save the info of the new minimum
                minAssocNum = refPointAssocNum[i];
                minAssocIndex = i;
            }
        }

        //Assign the next rarity index
        //First check to see if the reference point with the fewest associations has no associations
        if (refPointAssocNum[minAssocIndex] == 0) {
            //Set the reference points association number very high so it is not considered anymore
            refPointAssocNum[minAssocIndex] = INT_MAX;

        }        
        else {
            //Pick an adult from the points associated adults and assign the current rarity value
            //  Pick the final adult so it is easy to remove from the vector
            allAdults[refPointAssociations[minAssocIndex][refPointAssociations[minAssocIndex].size()-1]].rarity = rarityCount;

            //Remove the adult from the point's associations so it isn't assigned a new rarity
            refPointAssociations[minAssocIndex].pop_back();

            //Check to see if this was the last adult from this reference point
            if (refPointAssociations[minAssocIndex].size() == 0) {
                //No other adults are associated with the point, remove it from consideration
                refPointAssocNum[minAssocIndex] = INT_MAX;
            }
            else {
                //There are still adults so this point needs to be returned to
                //Add one to the points association count so other points can be considered
                refPointAssocNum[minAssocIndex]++; 
            }

            //Add one to the rarity counter for the next adult
            rarityCount++;
        }

        //Reset the trackers
        minAssocNum = allAdults.size() + 1;
        minAssocIndex = 0;
    }

    return;
}