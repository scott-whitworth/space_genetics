#define UNITTEST //should toggle unit testing functions throughout the code
#include "../Genetic_Algorithm/adult.h"
#include <iostream>
#include <vector>
#include <string> //allows us to use to_string and string
using std::cout;
using std::endl;

// Verifies that rankSort and rankDistanceSort are doing what we expect at a basic level (it can determine the better of two individuals)
// Compares two adults with set ranks and distances using rankSort and rankDistanceSort
// Uses the unitTesting constructor for Adult to make Adults with a known rank and distance so they can be better compared
// Prints cout statements to the terminal that tells you how the results of the comparison 
// Output: TRUE - passed all the tests as expected 
//         FALSE - on at least one test, at least one of the sorts did not behave as we expected
bool compareTwoAdults();

// Verifies that rankSort and rankDistanceSort will work on a vector of individuals
// Starts with two identical copies of a vector of adults and sorts each using rankSort or rankDistanceSort
// Prints cout statements to the terminal that show you what is going on with the sort
// Input: printStuff - TRUE -> will allow the user to either see all the individuals as they are organized 
//                     FALSE -> it will keep it from displaying all the cout statements
// Output: TRUE - All the vectors were sorted in the expected order
//         FALSE - During at least one of the tests, at least one of the vectors was not sorted how we expected
bool sortAdultVec(bool printStuff);

// Helper function for sortAdultVec to keep that function shorter
// Initializes a different vector of adults based on which test we are trying to perform
// Input: testNum - which test is being accomplished
//        a - the vector to be filled (should be empty)
//        print - whether or not it should print information about the vector to the terminal
// Output: all the current values are cleared from a and it is filled with adults corresponding to which test we are performing
bool differentTestsSetUp(int testNum, std::vector<Adult>& a, bool print);

// Ensures that all the tests got the expected results
// Holds the order in which the different Adults should be ordered and searches for issues with their actual order
// Input: testNum - the test we are performing
//        rs - Adults sorted using rankSort (should already be sorted when they're passed in)
//        rDS - Adults sorted using rankDistanceSort (should already be sorted when they're passed in)
// Output: returns TRUE if the vector was sorted correctly or FALSE if it was not
//         also prints cout statements explaining any discrepancies that were found
bool differentTestsCheck(int testNum, std::vector<Adult> & rs, std::vector<Adult>& rDS);

// Helper function for differentTestsCheck - cills correctRS and correctRDS with individuals in the order they should be depending on the test being executed
//      (based on testNum)
// Input: testNum - the test we are performing so it loads the correct vectors
//        correctRS & correctRDS - emptied and filled with the order the Adults should be in after each sort is performed
// Output: correctRS and correctRDS are filled with Adults in rankSort and rankDistanceSort order respectively for the test indicated by testNum
void loadCorrectOrders(int testNum, std::vector<Adult> & correctRS, std::vector<Adult> & correctRDS);

// Prints the rank, distance, and status of every element in the vector using unitTestingRankDistanceStatusPrint()
// Input: a - Adults to be printed 
// Output: prints the ranks, distances, and error statuses of each Adult in a
void vecPrinting(std::vector<Adult>& a);

#include "../Unit_Testing/testing_sorts.cpp"
