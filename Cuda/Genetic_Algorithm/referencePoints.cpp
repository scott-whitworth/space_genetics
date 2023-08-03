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
        //Add one to the total value to reflect the addition to the last index
        valTot += 1;

        //While loop to make sure the total value of the vector doesn't exceed the maximum value (based on cConstant's divisions)
        while (valTot > cConstants->divisions) {
            //Reset the value of the tracked index and subtract the value from the total value
            valTot -= newValues[indexTrack];
            newValues[indexTrack] = 0;

            //Go to the previous index and add one to its value
            indexTrack--;
            newValues[indexTrack] += 1;
            //add one to the total value
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

    //Add the last value the new point which is the number of divisions minus the total value, divided by the divisions
    //  This is because the total value of all divisions will equal divisions, so the last value is divisions minus the total value of the rest of the points
    //  It is then divided by divisions to make it a 0 to 1 scale 
    newPoint.push_back((cConstants->divisions - valTot)/cConstants->divisions);

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

            //If this is the initial call, the whole vector needs to be stored instead of calculating a determinant, so push the calculated determinant to the back of the norm vector
            if (initialMatrix) { 
                norm.push_back(pow(-1, i) * calcDet[0]); //only along the top row, j=0
            } 

            //Not the first call, calculate this part of the determinant
            else { 
                norm[0] += (pow(-1, i) * matrix[0][i] * calcDet[0]); 
            } 
        } 
    } 

    //Return the normal vector or the determinant
    return norm;
}

//Calculates the objective-by-objective relative cost for the adults passed into the function
void calculateRelCost (const cudaConstants *cConstants, std::mt19937_64 rng, ReferencePoints & refPoints, std::vector<Adult> & allAdults) {

    //Make sure the vectors are cleared
    refPoints.objWorst.clear();
    refPoints.objBest.clear();

    //Push back a random adult for each objective
    for (int obj = 0; obj < cConstants->missionObjectives.size(); obj++) {
        refPoints.objWorst.push_back(allAdults[(rng() % allAdults.size())]);
        refPoints.objBest.push_back(allAdults[(rng() % allAdults.size())]);
    }

    //Find the normalizations for each objective
    for (int obj = 0; obj < cConstants->missionObjectives.size(); obj++) {

        //Go through the rest of the adults to see if there are new best/worst values
        for (int indiv = 0; indiv < allAdults.size(); indiv++) {

            //Check to see if the adult has the new worst value for this objective
            //Because maximizations are stored as negatives, worse values are going to be larger values, regardless of objective
            if (allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) > refPoints.objWorst[obj].getParameters(cConstants->missionObjectives[obj])) {
                //Make sure this is not the same adult as the reference for finding the intercepts
                if (!duplicateCheck(allAdults[indiv], refPoints.objWorst[0], cConstants)) {
                    //Found a new worse value, store the new worst adult for the objective 
                    refPoints.objWorst[obj] = allAdults[indiv];
                }
            }
            //Check if the adult has the new best value for this objective that is not better than the dominationThreshold
            else if ((allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) < refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj])) 
                    && (allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) > cConstants->missionObjectives[obj].dominationThreshold)) {

                //Found a new best value, store the adult
                refPoints.objBest[obj] = allAdults[indiv];

                // std::cout<<"\n For obj "<<obj<<", new best by adult "<<indiv<<" of "<<allAdults[indiv].getParameters(cConstants->missionObjectives[obj])<<"\n";
            }          
        }
        // std::cout<<"\n For obj "<<obj<<", best adult has "<<refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj])<<"\n";
    }


    //Holds the points for the objective intercepts
    std::vector<double> intercepts;

    //Find intercepts if there are three objectives
    if (cConstants->missionObjectives.size() >= 2){
        //Get the on-plane vectors
        //Create a matrix of n, n-sized vectors between the worst points
        std::vector<std::vector<double>> matrix;

        //Calculate the on-plane vectors
        for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
            //i refers to the worst point on each objective
            //Add a new empty vector to the matrix
            matrix.push_back(std::vector<double>(0,0));

            for (int j = 0; j < cConstants->missionObjectives.size(); j++) {
                //j refers to each of the objectives of the ith point
                //The vectors are all relative to the first objective's worst adult's points
                //  so each component will be j adult's point minus the 0th adult's point
                //  the top row of the matrix will be all zeroes, but it won't matter for the calculation of the normal vector
                matrix[i].push_back(refPoints.objWorst[i].getParameters(cConstants->missionObjectives[j]) - refPoints.objWorst[0].getParameters(cConstants->missionObjectives[j]));
            }
        }

        //The matrix is created, calculate the plane's normal vector
        std::vector<double> normal = calcNormVector(matrix, true); 

        //Calculate the intercepts
        // The formula for intercept x (where the 0th adult is a, the normal objective is n, and the objectives are 1,2,3...) is (a1n1 + a2n2 + ...)/nx

        //Value to store the numerator for the intercepts
        double num = 0;

        //The numerator is consistent, use a for loop to calculate it
        for (int obj = 0; obj < cConstants->missionObjectives.size(); obj++) {
            //add the value of each component of the numerator to the total numerator value
            num += (refPoints.objWorst[0].getParameters(cConstants->missionObjectives[obj]) * normal[obj]);
        }

        //Calculate the intercepts
        for (int obj = 0; obj < cConstants->missionObjectives.size(); obj++) {
            //Divide the numerator by the intercept found for the objective
            intercepts.push_back(num / normal[obj]);

            //
            if(intercepts[obj] < refPoints.objWorst[obj].getParameters(cConstants->missionObjectives[obj])){
                intercepts[obj]=refPoints.objWorst[obj].getParameters(cConstants->missionObjectives[obj]);
            }
        }
    }

    //The normalization values have been found, go through the adults and calculate the objective costs
    for (int indiv = 0; indiv < allAdults.size(); indiv++) {

        //Calculate the cost for each objective
        for (int obj = 0; obj < cConstants->missionObjectives.size(); obj++){
            //Establish numerator and denominator for the normalization calculation
            //The numerator is the individual's objective value minus the best found objective value
            double translatedObj;

            // if((allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) - refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj])) > cConstants->missionObjectives[obj].equateTolerance){
            if(allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) > cConstants->missionObjectives[obj].dominationThreshold){
                translatedObj = allAdults[indiv].getParameters(cConstants->missionObjectives[obj]) - refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj]);
            }
            else{//If an individual's objective is lower than the domination threshold, then this is set to 0
                translatedObj = 0;
            }

            //The denominator should never be 0
            double denom=1;

            //If the normalization is based on objective intercepts (happens if there are two or more objectives), denominator is objective intercept - best found objective value
            if (cConstants->missionObjectives.size() >= 2){
                denom = intercepts[obj] - refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj]);
            }
            //If the normalization is not based on objective intercepts, the denominator is the worst objective value - the best objective value
            else {
                denom = refPoints.objWorst[obj].getParameters(cConstants->missionObjectives[obj]) - refPoints.objBest[obj].getParameters(cConstants->missionObjectives[obj]);
            }

            //Make sure there are no divide by 0 errors
            if (abs(denom) < cConstants->doublePrecThresh) {
            //if (denom < cConstants->missionObjectives[j].convergenceThreshold) {

                //minimization vs maximization handling
                if (cConstants->missionObjectives[obj].goal < 0) {
                    // denom = cConstants->missionObjectives[j].convergenceThreshold;
                    denom = cConstants->doublePrecThresh;
                }
                else {
                    // denom = -cConstants->missionObjectives[j].convergenceThreshold;
                    denom = -cConstants->doublePrecThresh;
                }
            }

            //Calculate the actual normailized objective
            allAdults[indiv].normalizedObj[obj] = translatedObj/denom;
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
void findAssociatedPoints (const ReferencePoints & refPoints, std::vector<Adult> & newAdults) {
    //Trackers for finding the closest reference point
    //In the normalized coordinates, the farthest corner (1,1,1) is at a distance of sqrt(3)=1.732... from the origin.
    double minDistVal = 2;      //Tracks the distance of the reference point currently the closest from the adult

    //Value that holds the distance from the adult to the reference point being checked
    double curDist = 2;

    //Find the closest point for all newAdults
    for (int i = 0; i < newAdults.size(); i++) {

        //For each adult check each reference point
        for (int j = 0; j < refPoints.points.size(); j++) {
            
            /* The following is to calculate the distance between the adult and the reference point itself (versus the origin line)

            //Calculate the distance to the reference point
            //First add the squares of the differences for each objective/component
            for (int k = 0; k < refPoints.points[0].size(); k++) {
                
                curDist += pow((refPoints.points[j][k] - newAdults[i].normalizedObj[k]), 2);
            }

            //Finally take the square root to find the distance
            curDist = sqrt(curDist);

            */

            //Calculate the distance from the adult to the line between the reference point and the origin
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
void calculateRarity (const ReferencePoints & refPoints, std::vector<Adult> & allAdults) {
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