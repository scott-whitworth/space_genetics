#include "../Genetic_Algorithm/referencePoints.h" //allows us access to the reference point functions

//gives each adult in the allAdults vector a rank
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

    //Used to store the current front of adults 
    //  first filled with the first front adults (best out of all population)
    //  filled with index of adults in allAdults
    std::vector<int> front;

    //2D vector stores which other adults each adult has dominated
    //1st dimension index for each adult
    //2nd dimension stores the indexes of other adults in allAdults that this adult has dominated
    std::vector<std::vector<int>> domination; 
    domination.resize(allAdults.size());

    //keep track of how many times an adult in oldAdults has been dominated by all other adult
    //each index corresponds to the same index within allAdults
    std::vector<int> dominatedByCount;
    
    //fill the vector with 0s to make sure the count is accurate
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through allAdults
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        //    becasue of mirroring aA[5] => aA[7] is the same as aA[7] => aA[5], only the upper half need to be checked
        for(int j = i+1; j < allAdults.size(); j++){
                        //Check to see if i dominates j
            if (dominationCheck(allAdults[i], allAdults[j], cConstants)){
                //Put the jth index in the set of individuals dominated by i
                domination[i].push_back(j);
                //add one to j's dominated by count
                dominatedByCount[j]++;
            }
            //Check to see if j dominates i
            else if (dominationCheck( allAdults[j], allAdults[i], cConstants)){
                //Put the ith index in the set of individuals dominated by j
                domination[j].push_back(i);
                //Add one to i's dominated by count
                dominatedByCount[i]++;
            }
            //implicit else: neither dominate eachother, don't do anything (don't update domination or dominatedByCount)
        }
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            front.push_back(i);
        }
    }

    //Used to assign rank number
    int rankNum = 1;
    //vector to store individuals' indexes in next front
    std::vector<int> newFront;

    //go until all individuals have been put in better ranks and none are left to give a ranking
    while(!front.empty()) {
        //empty the new front to put new individuals in
        std::vector<int>().swap(newFront);

        //loop through all individuals in old front
        //These individuals already have their rank set
        for(int k = 0; k < front.size(); k++){

            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < domination[front[k]].size(); l++){

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[front[k]][l]]--;
                
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[front[k]][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[front[k]][l]].rank = rankNum + 1;
                    //Add the index of the this adult to newFront
                    newFront.push_back(domination[front[k]][l]);                        
                }
            }
        }
        //increment the rank number
        rankNum++;

        //empty the current front
        std::vector<int>().swap(front);

        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
    }
}

//gives each adult in the allAdults a distance representing how different it is from other individuals
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants){
    
    //int that counts the number of the valid adults
    int validAdults = 0;

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(allAdults.begin(), allAdults.end(), rankSort);

    //checks if the adult is valid and then adds that index to the vector
    //the size of this vector will be used to find the distance for valid adults only
    for(int i = 0; i < allAdults.size(); i++){
        if(allAdults[i].errorStatus == VALID || allAdults[i].errorStatus == DUPLICATE){
            validAdults++;
        }
    }

    //set all to zero including the invalid adults
    for (int i = 0; i < allAdults.size(); i++ ){
        //reset each individual's distance
        allAdults[i].distance = 0.0;
    }

    //Iterate through the objectives
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
        
        //Sort allAdults based on the objective's target differences
        std::sort (allAdults.begin(), allAdults.begin()+validAdults, [i] (Adult a, Adult b) {return a.objTargetDiffs[i] < b.objTargetDiffs[i];});

        //Correct sort has been applied, apply the max distance to the extreme valid adults
        allAdults[0].distance += MAX_DISTANCE; 
        allAdults[validAdults-1].distance += MAX_DISTANCE;

        //For each individual besides the upper and lower bounds, make their distance equal to
        //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
        double normalParamLeft;
        double normalParamRight;

        //Tracks how many individuals have met the objective's threshold
        int metThreshold = 0;

        //Stores the value of the parameter which will be used to normalize the distances of the non-converged individuals
        //Value should be set to the largest valid value for the objective
        double normalizationValue = allAdults[validAdults-1].objTargetDiffs[i];

        // std::cout << "\nNormalization value: " << normalizationValue << "\n";
        
        //Add the distance to all non-converged individuals for this objective
        for(int j = metThreshold + 1; j < validAdults - 1; j++) {
            //Check to see if the adult has met the convergence threshold
            //Add the max distance to the adult if so
            if (allAdults[j].objTargetDiffs[i] < (cConstants->missionObjectives[i].allowedDifference + cConstants->missionObjectives[i].equateTolerance)) {
                allAdults[j].distance += MAX_DISTANCE;
            }
            //The adult has not met the convergence threshold for this objective, use the normal distance calculation
            else {
                //Divide left and right individuals by the largest value individual to normalize
                normalParamLeft = allAdults[j+1].objTargetDiffs[i] / normalizationValue;
                normalParamRight = allAdults[j-1].objTargetDiffs[i] / normalizationValue;

                // std::cout << "\tLeft: " << normalParamLeft << "\n\tRight: " << normalParamRight << "\n";

                //distance += abs((i+1) - (i-1))
                allAdults[j].distance += abs((normalParamLeft - normalParamRight));
            } 
        }
    }
}

//Assist function to determine the sort based on the selected algorithm
void mainSort(std::vector<Adult> & adults, const cudaConstants* cConstants, const int & sortSize) {
    //Switch statement determines the sort type based on the specified algorithm
    switch (static_cast<int>(cConstants->algorithm)) {
        //User wants a rank-rarity sort
        case RANK_RARITY:
            //sort by rank-rarity
            std::sort(adults.begin(), adults.begin()+sortSize, rankRaritySort);
            break;
        
        //Checking for unspecified as it will default to rank-distance
        case UNSPECIFIED:
            std::cout << "\nNo algorithm identified, defaulting to rank-distance!\n";

        //User wants rank-distance sorting
        case RANK_DISTANCE:
            //sort by rank distance
            std::sort(adults.begin(), adults.begin()+sortSize, rankDistanceSort);
            break;
    };   
}

//Not In Use
// Assist function for objectives that sorts an adult vector based on the goal 
void parameterSort(std::vector<Adult> & adults, const objective& sortObjective, const int & sortSize) {
    //Switch statement determines the type of sort needed to be done for the objective
    switch (static_cast<int>(sortObjective.goal)) {
        //Depending on the objective, the Adults must be sorted in a different order
        //While MIN_POS_DIFF and MIN_ORBITAL_POS_DIFF may seem quite similar, they are calcualted a little differently
        //  MIN_POS_DIFF is trying to make the position difference between the spacecraft and its target basically 0
        //  MIN_ORBITAL_POS_DIFF is trying to make the difference position difference between the spacecraft and its target equal to the orbital radius
        //  A similar thing is true with MIN_SPEED_DIFF & MIN_ORBITAL_SPEED_DIFF (going to 0 vs going to orbital speed)
        //  Because maximizations and minimizations are handled the same way, even maximizations need to be sorted by the lowest value
        case POS_DIFF:
            //Sort by lowest pos diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerPosDiff);
            break;

        case SPEED_DIFF:
            //Sort by lowest speed diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerSpeedDiff);
            break;

        case ORBIT_POS_DIFF:
            //Sort by lowest orbit pos diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerOrbitPosDiff);
            break;

        case ORBIT_SPEED_DIFF:
            //Sort by lowest orbit speed diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerOrbitSpeedDiff);
            break;

        case FUEL_SPENT:
            //Sort by the lowest fuel spent
            std::sort(adults.begin(), adults.begin()+sortSize, LowerFuelSpent);
            break;

        case TRIP_TIME:
            //sort by the lowest trip time
            std::sort(adults.begin(), adults.begin()+sortSize, LowerTripTime);
            break;

        case MARS_DIST:
            //Sort by lowest mars distance
            std::sort(adults.begin(), adults.begin()+sortSize, LowerMarsDist);
            break;

        case ORBIT_ASST:
            //Sort by lowest orbithChange
            std::sort(adults.begin(), adults.begin()+sortSize, LowerOrbitHChange);
            break;

        default:
            //No goal identified
            std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
            break;
    }
}