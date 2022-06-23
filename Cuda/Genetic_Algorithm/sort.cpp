//gives each adult in the allAdults vector a rank
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
        for(int j = 0; j < allAdults.size(); j++){
            //TODO: This should be j = i + 1
            //      This is will need to change all of the logic below VVV


            //TODO: Where is the best place to check for status conditions? In dominationCheck or here
            //check the status of both i and j and see if i is automatically dominated
            if(allAdults[i].errorStatus != VALID && allAdults[j].errorStatus == VALID){
                dominatedByCount[i]++;

            }//check the status of both i and j and see if j is automatically dominated
            else if(allAdults[j].errorStatus != VALID && allAdults[i].errorStatus == VALID){
                domination[i].push_back(j);
                
            }//if either both are valid or both are not valid, it will rank them normally
            //Check to see if i dominates j
            else if (dominationCheck(allAdults[i], allAdults[j], cConstants)){
                //Put the jth index in the set of individuals dominated by i
                //std::cout << "\n" << i << "th (i) Adult dominates " << j << "th (j) Adult!\n";
                domination[i].push_back(j);
            }
            //Check to see if j dominates i
            else if ( dominationCheck( allAdults[j], allAdults[i], cConstants) ){
                //std::cout << "\n" << j << "th (j) Adult dominates " << i << "th (i) Adult!\n";
                //Add one to i's dominated by count
                dominatedByCount[i]++; 
            }
            //Implicit else:
            //    (and i == j)
            //    If both are valid AND neither dominate eachother, then nothing changes (dominatedByCount / domination does not update)
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
    int validAdults;
    //starting rankSort to make sure nans are at the end of the array.
    std::sort(allAdults.begin(), allAdults.end(), rankSort);

    //checks if the adult is valid and then adds that index to the vector
    //the size of this vector will be used to find the distance for valid adults only
    for(int i = 0; i < allAdults.size(); i++){
        if(allAdults[i].errorStatus == VALID){
            validAdults++;
        }
    }

    //set all to zero including the invalid adults
    for (int i = 0; i < allAdults.size(); i++ ){
        //reset each individual's distance
        allAdults[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(allAdults.begin(), allAdults.begin() + validAdults, LowerPosDiff);
    //Set the boundaries
    allAdults[0].distance += MAX_DISTANCE; //+=1
    allAdults[validAdults - 1].distance += MAX_DISTANCE; //+=1


    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    for(int i = 1; i < validAdults - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffLeft = allAdults[i+1].posDiff/allAdults[validAdults - 1].posDiff;
        normalPosDiffRight = allAdults[i-1].posDiff/allAdults[validAdults - 1].posDiff;
        //distance = distance + abs((i+1) - (i-1))
        allAdults[i].distance = allAdults[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(allAdults[validAdults - 1].posDiff - allAdults[0].posDiff));
    }

    //Repeat above process for speedDiff
    if(cConstants->missionType == Rendezvous){//only do this for the rendezvous mission since it has 2 objectives
        std::sort(allAdults.begin(), allAdults.begin() + validAdults, LowerSpeedDiff);
        //Set the boundaries
        allAdults[0].distance += MAX_DISTANCE; //+=1
        allAdults[validAdults - 1].distance += MAX_DISTANCE; //+=1 //TODO:Was 0 for a bit, but isn't now consider which way is best

    
        //For each individual besides the upper and lower bounds, make their distance equal to
        //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
        double normalSpeedDiffLeft;
        double normalSpeedDiffRight;
        for(int i = 1; i < validAdults - 1; i++){
            //Divide left and right individuals by the worst individual to normalize
            normalSpeedDiffLeft = allAdults[i+1].speedDiff/allAdults[validAdults - 1].speedDiff;
            normalSpeedDiffRight = allAdults[i-1].speedDiff/allAdults[validAdults - 1].speedDiff;
            //distance = distance + abs((i+1) - (i-1))
            allAdults[i].distance = allAdults[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(allAdults[validAdults - 1].speedDiff - allAdults[0].speedDiff));
        }
    }
}
