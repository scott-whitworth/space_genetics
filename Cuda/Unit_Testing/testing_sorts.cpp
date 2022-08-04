bool compareTwoAdults(){
    //creates two adults to compare using the adult 
    Adult a1(1, 7);
    Adult a2(4, 9);

    //First comparison: Adult 1 (rank = 1 and distance = 7) & Adult 2 (rank = 4 and distance = 9)
    //compares two adults with different ranks and distances 
    if (!rankSort(a1, a2)){
        cout << "rankSort failed on the first comparison" << endl;
        return false;
    }
    if (!rankDistanceSort(a1, a2)){
        cout << "rankDistanceSort failed on the first comparison" << endl;
        return false;
    }

    //Second comparison: Adult 1 (rank = 1 and distance = 7) & Adult 3 (rank = 1 and distance = 6)
    //compares two adults with the same rank but different distances
    Adult a3(1,6);
    if (rankSort(a1, a3)){
        cout << "rankSort failed on the second comparison" << endl;
        return false;
    }
    if (!rankDistanceSort(a1, a3)){
        cout << "rankDistanceSort failed on the second comparison" << endl;
        return false;
    }

    //Third comparison: Adult 4 (rank = 1 and distance = 9) & Adult 1 (rank = 1 and distance = 7) where Adult 4 has an invalid status 
    //compares two adults with the same rank and different distances, but the one with the better distance is in an error state
    Adult a4(1, 9, OTHER_ERROR);
    if (rankSort(a4,a1)){
        cout << "rankSort failed on the third comparison" << endl;
        return false;
    }
    if (rankDistanceSort(a4,a1)){
        cout << "rankDistanceSort failed on the third comparison" << endl;
        return false;
    }

    //Fourth comparison: Adult 5 (rank = 2 and distance = 7) & Adult 1 (rank = 1 and distance = 7)   
    //compares two adults with the different ranks but the same distances 
    Adult a5(2, 7);
    if (rankSort(a5,a1)){
        cout << "rankSort failed om the fourth comparison" << endl;
        return false;
    }
    if (rankDistanceSort(a5,a1)){
        cout << "rankDistanceSort failed on the fourth comparison" << endl;
        return false;
    }

    //if all the tests passed, returns true
    return true;
}

bool sortAdultVec(bool printStuff){
    std::vector<Adult> forRankSort;
    std::vector<Adult> forRDSort; 
    
    //sortAdultVec tests sorting 6 different vectors of adults using rankSort and rankDistanceSort
    //these tests are set up in differentTestsSetUp and testNum allows it to start on the first test (test 1)
    int testNum = 1;

    //so far therse have not been any errors detected, so set this to true
    bool noErrors = true;

    //goes through this process for all the different vectors of Adults that are created by differentTestsSetUp
    //differentTestsSetUp will fill forRankSort with the correct values for each run (when a valid testNum is entered)
    //if the testNum is too large or too small it will return false and exit the while loop
    while (differentTestsSetUp(testNum, forRankSort, printStuff)) {
        forRDSort = forRankSort;
        
        if (printStuff){
            //Prints the adults before they are sorted -> just making sure forRDSort and forRankSort actually have same values ->this part can be deleted later
            cout << "Before Sorting (checking to make sure everything is correctly copied): " << endl;
            for (int j = 0; j < forRDSort.size(); j++){
                cout << forRDSort[j].unitTestingRankDistanceStatusPrint() << "; ";
            }
            cout << endl; 
        }

        //sorts the two vectors using their corresponding sorts
        std::sort(forRankSort.begin(), forRankSort.end(), rankSort);
        std::sort(forRDSort.begin(), forRDSort.end(), rankDistanceSort);

        if (printStuff){
            std::vector<Adult> goodRS;
            std::vector<Adult> goodRDS;
            loadCorrectOrders(testNum, goodRS, goodRDS);

            //Prints the results vs the correct version for rankSort
            cout << "rankSort (After Sorting): " << endl;
            vecPrinting(forRankSort);
            cout << "rankSort (Correct Version): " << endl;
            vecPrinting(goodRS);
            //Prints the results vs the correct version for rankDistanceSort
            cout << "rankDistanceSort (After Sorting): " << endl;
            vecPrinting(forRDSort);
            cout << "rankDistanceSort (Correct Version): " << endl;
            vecPrinting(goodRDS);
        }
        
        if (!differentTestsCheck(testNum,forRankSort,forRDSort)){
            noErrors = false;
        }

        //Scott had a problem with this because the check is pretty much the same as the rankSort function
        //Either delete this or leave it to "make sure getRank is working"
        //Checks to make sure that every element was organized correctly by rankSort
        for (int j = 0; j < forRankSort.size()-1; j++){
            if (forRankSort[j+1].errorStatus == VALID){ //if the next one is an error type, based on the way this is currently laid out, its rank might not actually be worse
                //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
                if (forRankSort[j].getRank() > forRankSort[j+1].getRank()){
                    cout << "There was a problem with rankSort on test " << testNum << endl;
                    noErrors = false;
                }
            }
        }
        //Checks to make sure every element was organized correctly by rankDistanceSort
        for (int k = 0; k < forRDSort.size()-1; k++){
            if (forRDSort[k+1].errorStatus == VALID){//if the next one is an error type, based on the way this is currently laid out, its rank might not actually be worse
                //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
                if (forRDSort[k].getRank() > forRDSort[k+1].getRank()){
                    cout << "There was a problem with rankDistanceSort on test " << testNum << endl;
                    noErrors = false;
                }
                //If the first and last have the same rank, but the smaller distance comes first, that is an issue
                else if (forRDSort[k].getRank() == forRDSort[k+1].getRank() && forRDSort[k].getDistance() < forRDSort[k+1].getDistance()){
                    cout << "There was a problem with rankDistanceSort on test " << testNum << endl;
                    noErrors = false;
                }
            }
        }
        testNum++; //increments testNum to prevent an infinite loop -> also moves on to the next test
    }
    //if it makes it through all the tests successfully, it passes and returns true
    return noErrors;
}

bool differentTestsSetUp(int testNum, std::vector<Adult>& a, bool print){
    if (testNum < 1 || testNum > 6){
        return false;
    }
    if (print){ 
        cout << "\nTest " << testNum << ": " << endl;
    }
    a.clear(); //starts by emptying the vector
    if(a.size() != 0){
        cout << "Clear did not work as expected " << endl;
    }

    //The first test is for a vector of 7 individuals each with a distinct rank and random distances
    if (testNum == 1){
        //Makes 7 Adults with ranks 1-7
        a.push_back(Adult(1,27, VALID));
        a.push_back(Adult(3,27, VALID));
        a.push_back(Adult(4,1, VALID));
        a.push_back(Adult(6,2,VALID));
        a.push_back(Adult(7,9));
        a.push_back(Adult(5,50));
        a.push_back(Adult(2,3));
    }
    //The second test looks at adults with tht same ranks, but differing distances
    else if (testNum == 2){
        //Makes 7 Adults that are all rank 1 with different distances 
        a.push_back(Adult(1,2, VALID));
        a.push_back(Adult(1,27, VALID));
        a.push_back(Adult(1,1, VALID));
        a.push_back(Adult(1,4, VALID));
        a.push_back(Adult(1,9, VALID));
        a.push_back(Adult(1,50, VALID));
        a.push_back(Adult(1,3, VALID));
    }
    //The third test looks at with different ranks but the same distance
    else if (testNum == 3){
        //Makes 7 Adults with different ranks but all the same distance (this should just be sorted by range) 
        a.push_back(Adult(8,11, VALID));
        a.push_back(Adult(3,11, VALID));
        a.push_back(Adult(11,11, VALID));
        a.push_back(Adult(2,11, VALID));
        a.push_back(Adult(6,11, VALID));
        a.push_back(Adult(5,11,VALID));
        a.push_back(Adult(1,11, VALID));
    }
    //The fourth test looks at adults who have a random ranks and distances (with some crossover in each rank and distance)
    else if (testNum == 4){
        //Makes 8 Adults with a variety of ranks and distances
        a.push_back(Adult(8,5, VALID));
        a.push_back(Adult(3,11, VALID));
        a.push_back(Adult(11,9, VALID));
        a.push_back(Adult(2,22, VALID));
        a.push_back(Adult(6,11, VALID));
        a.push_back(Adult(3,3, VALID));
        a.push_back(Adult(1,7, VALID));
        a.push_back(Adult(3,27, VALID));
    }
    //The fifth test examines duplicates and errors
    else if (testNum == 5){
        //Makes 5 Adults with three with all the same ranks and distances, but one is invalid and then two other values (one better and one worse)
        a.push_back(Adult(12,4, VALID));
        a.push_back(Adult(1,5, VALID));
        a.push_back(Adult(1,3, VALID));
        a.push_back(Adult(1,3, SUN_ERROR));
        a.push_back(Adult(1,3, VALID));
    }
    //The sixth test random ranks and distances with a combination of values with and without errors
    else if(testNum == 6){
        //Makes 10 Adults with a variety of ranks and distances, 4 have errors and 5 are valid
        a.push_back(Adult(8,5,VALID));
        a.push_back(Adult(3,11, SUN_ERROR));
        a.push_back(Adult(8,9,VALID));
        a.push_back(Adult(2,22,VALID));
        a.push_back(Adult(6,11, OTHER_ERROR));
        a.push_back(Adult(3,3,VALID));
        a.push_back(Adult(1,7, OTHER_ERROR));
        a.push_back(Adult(1,29, VALID));
        a.push_back(Adult(2,9, SUN_ERROR));
        a.push_back(Adult(3,3,VALID));
    }
    if (print){
        cout << "Before Sorting: ";
        vecPrinting(a);
    }

    return true;
}

//used to verify that the adults from each test were sorted correctly
bool differentTestsCheck(int testNum, std::vector<Adult> & rs, std::vector<Adult>& rDS){
    if (testNum < 1 || testNum > 6){ //if the testNum that gets passed in here cannot be evaluated, returns false
        cout << "testNum is invalid - exiting differentTestsCheck..." << endl;
        return false;
    }
    bool noErrorsInOrder = true;
    std::vector<Adult> correctRS;
    std::vector<Adult> correctRDS;
    loadCorrectOrders(testNum, correctRS, correctRDS); //loads vectors full adults in the order rankSort and rankDistanceSort should sort the adults
    //if the rank sorted and rank distance sorted vectors are not the same size, then there is a deeper issue
    //  with what is being passed in
    if (rs.size() != rDS.size()){
        cout << "A major issue occured and the vectors sorted using rankSort and rankDistanceSort are no longer the same lengths" << endl;
        return false;
    }

    if (rs.size() != correctRS.size() || rDS.size() != correctRDS.size()){
        cout << "There was an issue with loading one or more of the vectors" << endl;
        return false;
    }

    //now that we have confirmed the two vectors are the same length we can go through both vectors simulatenously
    for (int i = 0; i < rs.size(); i++){
        //confirms the correct individual is in each slot for the rank sorted individuals
        if (rs[i].rank != correctRS[i].rank || rs[i].distance != correctRS[i].distance || rs[i].errorStatus != correctRS[i].errorStatus){
            //prints a message if anything is in the wrong spot
            cout << "On test " << testNum << ", the vector sorted using rankSort has an incorrect value at index " << i << " of " << rs.size()-1;
            cout <<  ". \nThe values should be (" << correctRS[i].rank << "," << correctRS[i].distance << ",";
            if (correctRS[i].errorStatus == VALID){
                cout << "VALID) ";
            }
            else {
                cout << "invalid - error status) ";
            }
            cout << " but are (" << rs[i].rank << "," << rs[i].distance << ",";
            if (rs[i].errorStatus == VALID){
                cout << "VALID)." << endl;
            }
            else {
                cout << "invalid - error status)." << endl;
            }
            noErrorsInOrder = false;
        }
        //confirms each individual is in the correct spot for the rank distance sorted individuals
        if (rDS[i].rank != correctRDS[i].rank || rDS[i].distance != correctRDS[i].distance || rDS[i].errorStatus != correctRDS[i].errorStatus){
            //prints a message if anything is in the wrong spot
            cout << "On test " << testNum << ", the vector sorted using rankDistanceSort has an incorrect value at index " << i << " of " << rDS.size()-1;
            cout <<  ". \nThe values should be (" << correctRDS[i].rank << "," << correctRDS[i].distance << ",";
            if (correctRDS[i].errorStatus == VALID){
                cout << "VALID) ";
            }
            else {
                cout << "invalid - error status) ";
            }
            cout << " but are (" << rDS[i].rank << "," << rDS[i].distance << ",";
            if (rDS[i].errorStatus == VALID){
                cout << "VALID)." << endl;
            }
            else {
                cout << "invalid - error status)." << endl;
            }
            noErrorsInOrder = false;
        }
    }
    //if any errors were found this will return false, otherwise it will return true
    return noErrorsInOrder;
}

void loadCorrectOrders(int testNum, std::vector<Adult> & correctRS, std::vector<Adult> & correctRDS){
    correctRS.clear();
    correctRDS.clear();
    if (testNum < 1 || testNum > 6){
        cout << "testNum is invalid - cannot load the correct orders for the different vectors..." << endl;
        return;
    }
    else if (testNum == 1){
        //Fills correctRS with Adults in the order they should be after rankSort 
        //Should be sorted least to greatest rank
        correctRS.push_back(Adult(1,27, VALID));
        correctRS.push_back(Adult(2,3, VALID));
        correctRS.push_back(Adult(3,27, VALID));
        correctRS.push_back(Adult(4,1, VALID));
        correctRS.push_back(Adult(5,50, VALID));
        correctRS.push_back(Adult(6,2,VALID));
        correctRS.push_back(Adult(7,9, VALID));

        //Fills correctRDS with the Adults in the order they should be in after rankDistanceSort
        //because every Adult has a different rank, rankSort and rankDistanceSort should be in the same order
        correctRDS = correctRS;
    }
    else if (testNum == 2){
        //Fills correctRS with Adults in the order they should be in after rankSort
        //they all have the same rank and are all valid, so they should remain in the order they were entered in
        correctRS.push_back(Adult(1,2, VALID));
        correctRS.push_back(Adult(1,27, VALID));
        correctRS.push_back(Adult(1,1, VALID));
        correctRS.push_back(Adult(1,4, VALID));
        correctRS.push_back(Adult(1,9, VALID));
        correctRS.push_back(Adult(1,50, VALID));
        correctRS.push_back(Adult(1,3, VALID));

        //Fills correctRDS with Adults in the order they should be in after rankDistanceSort
        //because every individual has the same rank, their order is soley dependent on their distances (greatest to least)
        correctRDS.push_back(Adult(1,50, VALID));
        correctRDS.push_back(Adult(1,27, VALID));
        correctRDS.push_back(Adult(1,9, VALID));
        correctRDS.push_back(Adult(1,4, VALID));
        correctRDS.push_back(Adult(1,3, VALID));
        correctRDS.push_back(Adult(1,2, VALID));
        correctRDS.push_back(Adult(1,1, VALID));
    }
    else if (testNum == 3){
        //Fills correctRS with Adults in the order they should be in after rankSort
        //Much like test 1, these should only be sorted based on rank because there are no repeated ranks and all the ranks vary 
        //(only distances are the same)
        correctRS.push_back(Adult(1,11, VALID));
        correctRS.push_back(Adult(2,11, VALID));
        correctRS.push_back(Adult(3,11, VALID));
        correctRS.push_back(Adult(5,11,VALID));
        correctRS.push_back(Adult(6,11, VALID));
        correctRS.push_back(Adult(8,11, VALID));
        correctRS.push_back(Adult(11,11, VALID));

        //Fills correctRDS with Adults in the order they should be in after rankDistanceSort
        //Much like test 1, these should only be sorted based on rank because there are no repeated ranks and all the ranks vary 
        //Thus, the order of correctRS and correctRDS should be the same
        correctRDS.push_back(Adult(1,11, VALID));
        correctRDS.push_back(Adult(2,11, VALID));
        correctRDS.push_back(Adult(3,11, VALID));
        correctRDS.push_back(Adult(5,11,VALID));
        correctRDS.push_back(Adult(6,11, VALID));
        correctRDS.push_back(Adult(8,11, VALID));
        correctRDS.push_back(Adult(11,11, VALID));
    }
    else if (testNum == 4){
        //Fills correctRS with Adults in the order they should be in after rankSort
        //There is a variety of ranks and distances, so Adults who have the same rank will appear in the order they were entered
        correctRS.push_back(Adult(1,7, VALID));
        correctRS.push_back(Adult(2,22, VALID));
        correctRS.push_back(Adult(3,11, VALID));
        correctRS.push_back(Adult(3,3, VALID));
        correctRS.push_back(Adult(3,27, VALID));
        correctRS.push_back(Adult(6,11, VALID));
        correctRS.push_back(Adult(8,5, VALID));
        correctRS.push_back(Adult(11,9, VALID));

        //Fills correctRDS with Adults in the order they should be in after rankDistanceSort
        //This will be similar to the order for rankSort, but not the same because Adults that share a rank will also be sorted by distance
        correctRDS.push_back(Adult(1,7, VALID));
        correctRDS.push_back(Adult(2,22, VALID));
        correctRDS.push_back(Adult(3,27, VALID));
        correctRDS.push_back(Adult(3,11, VALID));
        correctRDS.push_back(Adult(3,3, VALID));
        correctRDS.push_back(Adult(6,11, VALID));
        correctRDS.push_back(Adult(8,5, VALID));
        correctRDS.push_back(Adult(11,9, VALID));
    }
    else if (testNum == 5){ //Test 5 deals with copies of an Adult, including one whose status is set to invalid
        //Fills correctRS with Adults in the order they should be in after rankSort
        //The individual with an invalid status should be put at the end of the vector, and the others will be sorted solely based on rank
        correctRS.push_back(Adult(1,5, VALID));
        correctRS.push_back(Adult(1,3, VALID));
        correctRS.push_back(Adult(1,3, VALID));
        correctRS.push_back(Adult(12,4, VALID));
        correctRS.push_back(Adult(1,3, SUN_ERROR));

        //Fills correctRDS with Adults in the order they should be in after rankDistanceSort
        //Thanks to the order the Adults were entered into the vector in, rankDistanceSort and rankSort happen to have the same values
        //If one of the Adults with a VALID errorStatus and rank 1 distance 3 was entered before the rank 1 distance 5 individual
        //then correctRDS would be different than correctRS
        correctRDS = correctRS;
    }
    else if (testNum == 6){ //10 Adults with a variety of ranks, distances, and errorStatuses (also some duplicates)
        //Fills correctRS with Adults in the order they should be in after rankSort
        //The Adults should be ordered based on their validity, then their rank
        //Those that are invalid are not sorted by rank, but appear in the order they were entered in
        correctRS.push_back(Adult(1,29, VALID));
        correctRS.push_back(Adult(2,22,VALID));
        correctRS.push_back(Adult(3,3,VALID));
        correctRS.push_back(Adult(3,3,VALID));
        correctRS.push_back(Adult(8,5,VALID));
        correctRS.push_back(Adult(8,9,VALID));
        correctRS.push_back(Adult(3,11, SUN_ERROR));
        correctRS.push_back(Adult(6,11, OTHER_ERROR));
        correctRS.push_back(Adult(1,7, OTHER_ERROR));
        correctRS.push_back(Adult(2,9, SUN_ERROR));

        //Fills correctRDS with Adults in the order they should be in after rankDistanceSort
        //The Adults should be ordered based on their validity, then their rank, then their distance (greatest to least)
        //Those that are invalid are not sorted by rank or distance, but appear in the order they were entered in
        correctRDS.push_back(Adult(1,29, VALID));
        correctRDS.push_back(Adult(2,22,VALID));
        correctRDS.push_back(Adult(3,3,VALID));
        correctRDS.push_back(Adult(3,3,VALID));
        correctRDS.push_back(Adult(8,9,VALID));
        correctRDS.push_back(Adult(8,5,VALID));
        correctRDS.push_back(Adult(3,11, SUN_ERROR));
        correctRDS.push_back(Adult(6,11, OTHER_ERROR));
        correctRDS.push_back(Adult(1,7, OTHER_ERROR));
        correctRDS.push_back(Adult(2,9, SUN_ERROR));
    }
}

//prints the rank, distance, and errorStatus of each element of the vector
void vecPrinting(std::vector<Adult>& a){
    for (int i = 0; i < a.size(); i++){
        cout << a[i].unitTestingRankDistanceStatusPrint() << "; ";
    }
    cout << endl;
}

