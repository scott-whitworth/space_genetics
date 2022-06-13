#define UNITTEST //should toggle unit testing functions throughout the code
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
void vecPrinting(std::vector<Adult>& a);

#include "../Unit_Testing/testing_sorts.cpp"