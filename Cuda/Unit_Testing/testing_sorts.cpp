#include "../Genetic_Algorithm/adult.h"
#include <iostream>
#include <vector>
#include <string> //allows us to use to_string and string
using std::cout;
using std::endl;

// This function compares two adults with set ranks and distances using rankSort and rankDistanceSort
// Uses the unitTesting constructor for Adult to make Adults with a known rank and distance so they can be better compared
// Prints cout statements to the terminal that tells you how the results of the comparison 
// Output: Returns true if it passed all the tests as expected and exits and returns false if something went wrong
bool compareTwoAdults();

// This function uses a combination of std::sort with rankSort and rankDistance sort
// Starts with two identical copies of a vector of adults and sorts them using rankSort and rankDistanceSort
// Prints cout statements to the terminal that show you what is going on with the sort
// Returns a bool that says whether or not it passed all the tests
bool sortAdultVec();

// This function is just here to keep sortAdultVec from being to insanely long - it initializes a different vector of adults based on which test we are trying to perform
// Input: which test is being accomplished
// Output: all the current values are cleared from a and it is filled with adults corresponding to which test we are performing
void differentTestsSetUp(int testNum, std::vector<Adult>& a);

//prints the rank, distance, and status of every element in the vector using unitTestingRankDistanceStatusPrint()
// Input: a vector of adults
// Output: prints the ranks, distances, and error statuses of each adult
void vecPrinting(std::vector<Adult>& a);

// this function will check if giveRank is functioning as expected
bool dominationTest();

void giveRankTest(std::vector<Adult> & allAdults);

int test_main(const cudaConstants* cConstants){
    if (compareTwoAdults()){ //if it sucessfully makes it through compareTwoAdults, it prints a confirmation message for the user 
        cout << "PASSED: differentiated the better of two adults" << endl;
    }
    if (sortAdultVec()){ //if it successfully makes it through sortAdultVec, it tells the user it was sucessful
      cout << "PASSED: sorted adults as expected" << endl;
    }
    if(dominationTest()){
        cout << "dominationTest has ended" <<endl;
    }
    return 0;
}

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

bool sortAdultVec(){
    std::vector<Adult> forRankSort;
    std::vector<Adult> forRDSort; 
    
    int testTypes = 6; //a number to represent how many different tests you want to perform (should correspond to the number of available options in differentTestSetUps)
    
    //goes through this process for all the different vectors of Adults that are created by differentTestsSetUp
    for (int i = 1; i <= testTypes; i++) {
        cout << "\n Test " << i << ": " << endl;
        differentTestsSetUp(i, forRankSort); //depending on which time through the loop this is, sets up fills forRankSort with values 
        forRDSort = forRankSort;
        
        //Prints the adults before they are sorted -> just making sure forRDSort and forRankSort actually have same values ->this part can be deleted later
        cout << "Before Sorting: ";
        for (int j = 0; j < forRDSort.size(); j++){
            cout << forRDSort[j].unitTestingRankDistanceStatusPrint() << "; ";
        }
        cout << endl;

        //sorts the two vectors using their corresponding sorts
        std::sort(forRankSort.begin(), forRankSort.end(), rankSort);
        std::sort(forRDSort.begin(), forRDSort.end(), rankDistanceSort);
        
        //Prints the results of the sorts
        cout << "After Sorting (rankSort): " << endl;
        vecPrinting(forRankSort);
        cout << "After Sorting (rankDistanceSort): " << endl;
        vecPrinting(forRDSort);
        
        //Once I take out most of the cout statements, have it go through every element, comparing neighbors to ensure their ranks and distances have the right relationship
        //if the last one is an error tyoe, based on the way this is currently laid out, its rank might not actually be worse
        if (forRankSort[forRankSort.size()-1].errorStatus == VALID){
            //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
            if (forRankSort[0].getRank() > forRankSort[forRankSort.size()-1].getRank()){
                cout << "There was a problem with rankSort on test " << i << endl;
                return false;
            }
        }
        //if the last one is an error tyoe, based on the way this is currently laid out, its rank might not actually be worse
        if (forRDSort[forRDSort.size()-1].errorStatus == VALID){
            //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
            if (forRDSort[0].getRank() > forRDSort[forRDSort.size()-1].getRank()){
                cout << "There was a problem with rankDistanceSort on test " << i << endl;
                return false;
            }
            //If the first and last have the same rank, but the smaller distance comes first, that is an issue
            else if (forRDSort[0].getRank() == forRDSort[forRDSort.size()-1].getRank() && forRDSort[0].getDistance() < forRDSort[forRDSort.size()-1].getDistance()){
                cout << "There was a problem with rankDistanceSort on test " << i << endl;
                return false;
            }
        }
    }
    return true;
}

void differentTestsSetUp(int testNum, std::vector<Adult>& a){
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
        a.push_back(Adult(7,9,VALID));
        a.push_back(Adult(5,50,VALID));
        a.push_back(Adult(2,3,VALID));
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
    //The fourth test looks at adults who have a randp, ranks and distances (with some crossover in each rank and distance)
    else if (testNum == 4){
        //Makes 7 Adults with a variety of ranks and distances
        a.push_back(Adult(8,5, VALID));
        a.push_back(Adult(3,11, VALID));
        a.push_back(Adult(11,9, VALID));
        a.push_back(Adult(2,22, VALID));
        a.push_back(Adult(6,11, VALID));
        a.push_back(Adult(3,3, VALID));
        a.push_back(Adult(1,7, VALID));
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
        //Makes 9 Adults with a variety of ranks and distances, 4 have errors and 5 are valid
        a.push_back(Adult(8,5,VALID));
        a.push_back(Adult(3,11, SUN_ERROR));
        a.push_back(Adult(8,9,VALID));
        a.push_back(Adult(2,22,VALID));
        a.push_back(Adult(6,11, OTHER_ERROR));
        a.push_back(Adult(3,3,VALID));
        a.push_back(Adult(1,7, OTHER_ERROR));
        a.push_back(Adult(1,29, VALID));
        a.push_back(Adult(2,9, SUN_ERROR));
    }
    cout << "Before Sorting: ";
    vecPrinting(a);
}

//prints the rank, distance, and errorStatus of each element of the vector
void vecPrinting(std::vector<Adult>& a){
    for (int i = 0; i < a.size(); i++){
        cout << a[i].unitTestingRankDistanceStatusPrint() << "; ";
    }
    cout << endl;
}

bool dominationTest(){

    //vector for the test
    std::vector<Adult> GRtest;

    int r = 1;//rank
    int d = 1;//distance
    int vectSize = 6;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        GRtest.push_back(Adult(r, d)); 
    }
    
    //giving each adult specific posDiffs and speedDiffs for expected results
    //should be rank 5; expected dominatedByCount = 4
    GRtest[0].posDiff = 14;
    GRtest[0].speedDiff = 14;
    //should be rank 4; expected dominatedByCount = 2
    GRtest[1].posDiff = 10;
    GRtest[1].speedDiff = 10;
    //should be rank 3; expected dominatedByCount = 1?
    GRtest[2].posDiff = 8;
    GRtest[2].speedDiff = 8;
    //should be rank 2; expected dominatedByCount = 1
    GRtest[3].posDiff = 7;
    GRtest[3].speedDiff = 7;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[4].posDiff = 3;
    GRtest[4].speedDiff = 3;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[5].posDiff = 1;
    GRtest[5].speedDiff = 1;
    
    giveRankTest(GRtest);
    cout << "Ranks:" << endl;
    for(int i = 0; i < GRtest.size(); i++){
        cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << endl;  
    }
    
    return true;
    
}
void giveRankTest(std::vector<Adult> & allAdults) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf
    
    //Used to store the current front of adults. first filled with the first front adults(best out of all population)
    // filled with index of adults in allAdults
    std::vector<int> front;

    //TODO: Pull Adult::dominates and Adult::dominatedByCount into this function
    // probably a 2D vector, or an array of vectors
    //This 2D vector will store which other adults each adult has dominated
    //1st dimension will be a spot for each adult in allAdults
    //2nd dimension will store the indexes of adults in allAdults that the adult in the 1st dimension has dominated
    std::vector<std::vector<int>> domination; 

    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    //Note: fill the vector with 0s to make sure the count is accurate
    //TODO: unit test to make sure the whole vector is actually initially filled with 0's and not just the first index or the original vector size
    std::vector<int> dominatedByCount(0); 

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){
        cout << "i: " << i << " ";
        
        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){
            cout << " j: " << j << " ";

            //Check to see if i dominates j
            if (dominationCheckTest(allAdults[i], allAdults[j])){
                //Put the jth index in the set of individuals dominated by i
                domination[i].push_back(j);
                cout << " i dominated j ";
                //TODO: will this add too many to j's domination count? When it i's current value reaches j's current value it will have already recorded the dominaton here, but it will be recorded again 
                //Add one to j's dominated by count
                //dominatedByCount[j]++; 
            }
            //Check to see if j dominates i
            else if (dominationCheckTest(allAdults[j], allAdults[i])) {
                //TODO: this may have the same redundancy that was mentioned above with things being added to this vector too many times
                //Put the ith index in the set of individuals dominated by j
                //domination[j].push_back(i);
                cout << " j dominated i ";
                //Add one to i's dominated by count
                dominatedByCount[i]++; 
            }
            
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            front.push_back(i);
        }
        cout << "i was dominated by j " << dominatedByCount[i] << " times ";
        cout << endl;
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
            for(int l = 0; l < domination[k].size(); l++){

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[k][l]]--;
                cout << "index: " << domination[k][l] << " count: " << dominatedByCount[domination[k][l]] << " ";
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[k][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[k][l]].rank = rankNum + 1;
                    cout << "index: " << domination[k][l] << " current rank: " << allAdults[domination[k][l]].rank << " ";
                    //Add the index of the this adult to newFront
                    newFront.push_back(l);                        
                }
                cout << endl;
            }
        }
        //increment the rank number
        rankNum++;
        cout << " current rank: " << rankNum << endl;
        //empty the current front
        std::vector<int>().swap(front);

        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
    }
}
