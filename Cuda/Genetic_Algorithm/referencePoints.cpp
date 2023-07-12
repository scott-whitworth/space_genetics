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

//Calculates the objective-by-objective relative cost for the adults passed into the function
void calculateRelCost (const cudaConstants *cConstants, std::vector<Adult> & allAdults) {
    //Vector will hold the relative normalization values for all of the objectives
    //  The normalization value should be the overall worst value for the objective
    std::vector<double> normalizations (cConstants->missionObjectives.size(), 0);

    //Find the normalizations for each objective
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //First set the normalization value to first adult in the allAdult vector
        normalizations[i] = allAdults[0].getParameters(cConstants->missionObjectives[i]);

        //Go through the rest of the adults to see if there is a worst value
        for (int j = 1; j < allAdults.size(); j++){

            //Because maximizations are stored as negatives, worse values are going to be larger values, regardless of objective
            //Check to see if this adult has a worse value
            if (allAdults[j].getParameters(cConstants->missionObjectives[i]) > normalizations[i]) {

                //Found a worse value, store it as a new normalization
                normalizations[i] = allAdults[j].getParameters(cConstants->missionObjectives[i]);
            }
        }
    }

    //The normalization values have been found, go through the adults and calculate the objective costs
    for (int i = 0; i < allAdults.size(); i++) {

        //Calculate the cost for each objective
        for (int j = 0; j < cConstants->missionObjectives.size(); j++){
            //Store this objective's progress with the calculation (adult's objective value / normalization)
            double objProg = (allAdults[i].getParameters(cConstants->missionObjectives[j])/normalizations[j]);

            //If the objective is a maximization, the progress is not a 0 to 1 scale, so we need the inverse
            if (cConstants->missionObjectives[j].goal > 0) {
                objProg = 1/objProg;
            }

            //Calculate the cost as 1 - objProg
            allAdults[i].objectiveCost[j] = 1-objProg;
        }
    }

    return;
}

//Will find the closest reference points to each adult in the newAdults vector
void findAssociatedPoints (const cudaConstants *cConstants, const ReferencePoints & refPoints, std::vector<Adult> & newAdults) {
    //Trackers for finding the closest reference point
    double minDistVal = 2;      //Tracks the distance of the reference point currently the closest from the adult

    //Value that holds the distance from the adult to the reference point being checked
    double curDist = 2;

    //Find the closest point for all newAdults
    for (int i = 0; i < newAdults.size(); i++) {

        //For each adult check each reference point
        for (int j = 0; j < refPoints.points.size(); j++) {
            
            //Calculate the distance to the reference point
            //First add the squares of the differences for each objective/component
            for (int k = 0; k < refPoints.points[0].size(); k++) {
                
                curDist += pow((refPoints.points[j][k] - newAdults[i].objectiveCost[k]), 2);
            }

            //Finally take the square root to find the distance
            curDist = sqrt(curDist);

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
                refPointAssocNum[minAssocIndex]--; 
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