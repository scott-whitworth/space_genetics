//gives each adult in the allAdults vector a rank
//Branching test.
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

    //Used to store the current front of adults. first filled with the first front adults(best out of all population)
    // filled with index of adults in allAdults
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

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        for(int j = i+1; j < allAdults.size(); j++){

            //TODO: Where is the best place to check for status conditions? In dominationCheck or here
            //check the status of both i and j and see if i is automatically dominated
            /*
            if(allAdults[i].errorStatus != VALID && allAdults[j].errorStatus == VALID){
                dominatedByCount[i]++;
                domination[j].push_back(i);
            }//check the status of both i and j and see if j is automatically dominated
            else if(allAdults[j].errorStatus != VALID && allAdults[i].errorStatus == VALID){
                domination[i].push_back(j);
                dominatedByCount[j]++;
                
            }
            */
            //if either both are valid or both are not valid, it will rank them normally
            //Check to see if i dominates j
            if (dominationCheck(allAdults[i], allAdults[j], cConstants)){
                //Put the jth index in the set of individuals dominated by i
                //std::cout << "\n" << i << "th (i) Adult dominates " << j << "th (j) Adult!\n";
                domination[i].push_back(j);
                dominatedByCount[j]++;
            }
            //Check to see if j dominates i
            else if ( dominationCheck( allAdults[j], allAdults[i], cConstants) ){
                //std::cout << "\n" << j << "th (j) Adult dominates " << i << "th (i) Adult!\n";
                //Add one to i's dominated by count
                dominatedByCount[i]++;
                domination[j].push_back(i);
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            front.push_back(i);

            //std::cout << "\n\nAdult #" << i << " ranked " << 1;
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

                    //std::cout << "\n\nAdult #" << domination[front[k]][l] << " ranked " << rankNum+1;

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
    //std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~FINISHED RANKING~~~~~~~~~~~~~~~~~~~~~~~~\n";
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
        
        //Sort allAdults based on the objective's goal
        parameterSort(allAdults, cConstants->missionObjectives[i], validAdults); 

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
        //Value should be set to the largest valid value for the objective, which will be located at different points depending on the objective
        double normalizationValue = 0;

        //Find how many adults have met the threshold
        //Need to determine the optimization direction for the right comparison operator
        if (cConstants->missionObjectives[i].goal < 0) { //Minimize
            //For minimizations, the normalization value will be set to the last valid value
            normalizationValue = allAdults[validAdults-1].getParameters(cConstants->missionObjectives[i]);

            //Go through the non-extreme valid adults and see if they have met the threshold
            for (int j = 1; j < validAdults-1; j++) {
                if (allAdults[j].getParameters(cConstants->missionObjectives[i]) < cConstants->missionObjectives[i].convergenceThreshold) {
                    //Add one to the metThreshold tracker
                    metThreshold++;

                    //Add the max distance to the adult's distance to further promote it
                    allAdults[j].distance += MAX_DISTANCE;
                }
                else {
                    //Since the allAdults vector is sorted for this objective, if the code is here, it means all adults who meet the threshold have been passed
                    //So it is safe to break the for loop
                    break;
                }
            }
            
        }
        else if (cConstants->missionObjectives[i].goal > 0) { //Maximize
            //Go through the non-extreme valid adults and see if they have met the threshold
            for (int j = 1; j < validAdults-1; j++) {
                if (allAdults[j].getParameters(cConstants->missionObjectives[i]) > cConstants->missionObjectives[i].convergenceThreshold) {
                    //Add one to the metThreshold tracker
                    metThreshold++;

                    //Add the max distance to the adult's distance to further promote it
                    allAdults[j].distance += MAX_DISTANCE;
                }
                else {
                    //Since the allAdults vector is sorted for this objective, if the code is here, it means all adults who meet the threshold have been passed
                    //First, set the normalization value to the highest non-converged value, which is the metThreshold'th individual
                    normalizationValue = allAdults[metThreshold].getParameters(cConstants->missionObjectives[i]);

                    //Next it is safe to break the for loop
                    break;
                }
            }
        }
        else {
            //Error getting the goal
            std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        }
        
        //Add the distance to all non-converged individuals for this objective
        for(int j = metThreshold + 1; j < validAdults - 1; j++) {
            //Divide left and right individuals by the largest value individual to normalize
            normalParamLeft = allAdults[j+1].getParameters(cConstants->missionObjectives[i]) / normalizationValue;
            normalParamRight = allAdults[j-1].getParameters(cConstants->missionObjectives[i]) / normalizationValue;

            //distance += abs((i+1) - (i-1))
            allAdults[j].distance += abs((normalParamLeft - normalParamRight));
        }
    }
}

// Assist function for objectives that sorts an adult vector based on the goal 
void parameterSort(std::vector<Adult> & adults, const objective& sortObjective, const int & sortSize) {
    //Switch statement determines the type of sort needed to be done for the objective
    switch (static_cast<int>(sortObjective.goal)) {
        case MIN_POS_DIFF:
            //Sort by lowest pos diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerPosDiff);
            break;

        case MIN_SPEED_DIFF:
            //Sort by lowest speed diff
            std::sort(adults.begin(), adults.begin()+sortSize, LowerSpeedDiff);
            break;

        case MIN_FUEL_SPENT:
            //Sort by the lowest fuel spent
            std::sort(adults.begin(), adults.begin()+sortSize, LowerFuelSpent);
            break;

        case MIN_TRIP_TIME:
            //sort by the lowest trip time
            std::sort(adults.begin(), adults.begin()+sortSize, LowerTripTime);
            break;

        case MAX_SPEED_DIFF:
            //Sort by maximum speed diff
            std::sort(adults.begin(), adults.begin()+sortSize, HigherSpeedDiff);
            break;

        default:
            //No goal identified
            std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
            break;
    }
}