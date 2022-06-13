/*#define UNITTEST //should toggle unit testing functions throughout the code
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
// Input: the bool printStuff will allow the user to either see all the individuals as they are organized (true) or it will keep it from displaying all the cout statements (false)
// Returns a bool that says whether or not it passed all the tests
bool sortAdultVec(bool printStuff);

// This function is just here to keep sortAdultVec from being to insanely long - it initializes a different vector of adults based on which test we are trying to perform
// Input: which test is being accomplished, the vector to be filled, and whether or not it should print information about the vector to the terminal
// Output: all the current values are cleared from a and it is filled with adults corresponding to which test we are performing
void differentTestsSetUp(int testNum, std::vector<Adult>& a, bool print);

//prints the rank, distance, and status of every element in the vector using unitTestingRankDistanceStatusPrint()
// Input: a vector of adults
// Output: prints the ranks, distances, and error statuses of each adult
void vecPrinting(std::vector<Adult>& a);*/

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
    
    int testTypes = 6; //a number to represent how many different tests you want to perform (should correspond to the number of available options in differentTestSetUps)
    //TODO: This is a little weird to set this here, just because it is easy to miss / forget to change
    //      I might suggest either having differentTestSetUp return false if the entered index is not valid 
    //          That way you can have a while(differentTestsSetUp) loop without needing to 'export' the number of test types
    //      or make testTypes a const global var inside testing_rank (right above the function prototype) TEST_TYPE_NUM
    //          I don't love this, just because it is a kind of unnecessary variable given the above pattern
    //
    // Also, this comment could use some work:
    // int testTypes = 6; // Total number of tests to be performed (should correspond to actual options in differentTestsSetUp)
    // to pull apart what is there: 'a number to represent' is kind of obvious, all int variables are numbers that represent something
    //                              long comments are hard to read, you can either 1) make them more concise
    //                                                                          or 2) use more than one line
    
    //goes through this process for all the different vectors of Adults that are created by differentTestsSetUp
    for (int i = 1; i <= testTypes; i++) {
        if (printStuff){ 
            cout << "\nTest " << i << ": " << endl;
        }
        differentTestsSetUp(i, forRankSort, printStuff); //depending on which time through the loop this is, sets up fills forRankSort with values 
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
            //Prints the results of the sorts
            cout << "After Sorting (rankSort): " << endl;
            vecPrinting(forRankSort);
            cout << "After Sorting (rankDistanceSort): " << endl;
            vecPrinting(forRDSort);
        }
        
        //Checks to make sure that every element was organized correctly by rankSort
        for (int j = 0; j < forRankSort.size()-1; j++){
            if (forRankSort[j+1].errorStatus == VALID){ //if the next one is an error type, based on the way this is currently laid out, its rank might not actually be worse
                //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
                //TODO: There is a bit of a tautology going on here, rankSort is functionally using getRank
                //      This is a hard one to get around as it is such a basic part of the class
                //      I might suggest having a verified order (like hand calculate what the order should be), then check against that
                if (forRankSort[j].getRank() > forRankSort[j+1].getRank()){
                    cout << "There was a problem with rankSort on test " << i << endl;
                    return false;
                }
            }
        }

        //TODO: Same issue as outlined above
        //      You are using the same methods to confirm this is working as the methods used in the thing you are testing
        //      A different way to address this would be to methodically check the underlying methods before running this code (which you might already be doing)

        //Checks to make sure every element was organized correctly by rankDistanceSort
        for (int k = 0; k < forRDSort.size()-1; k++){
            if (forRDSort[k+1].errorStatus == VALID){//if the next one is an error type, based on the way this is currently laid out, its rank might not actually be worse
                //If the rank of the last thing is lower/better than that of the first, there is definitely something going wrong
                if (forRDSort[k].getRank() > forRDSort[k+1].getRank()){
                    cout << "There was a problem with rankDistanceSort on test " << i << endl;
                    return false;
                }
                //If the first and last have the same rank, but the smaller distance comes first, that is an issue
                else if (forRDSort[k].getRank() == forRDSort[k+1].getRank() && forRDSort[k].getDistance() < forRDSort[k+1].getDistance()){
                    cout << "There was a problem with rankDistanceSort on test " << i << endl;
                    return false;
                }
            }
        }
    }
    //if it makes it through all the tests successfully, it passes and returns true
    return true;
}

void differentTestsSetUp(int testNum, std::vector<Adult>& a, bool print){
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
    if (print){
        cout << "Before Sorting: ";
        vecPrinting(a);
    }
}

//prints the rank, distance, and errorStatus of each element of the vector
void vecPrinting(std::vector<Adult>& a){
    for (int i = 0; i < a.size(); i++){
        cout << a[i].unitTestingRankDistanceStatusPrint() << "; ";
    }
    cout << endl;
}

