#define UNITTEST
#include <iostream>
using std::cout;
using std::endl;
#include "../Unit_Testing/testing_sorts.h"
#include "../Unit_Testing/testing_genetics.h"

int main(){
    bool lotsOfCouts = true; //determines whether or not all adults' ranks, distances, and statuses will be printed
    if (compareTwoAdults()){ //if it sucessfully makes it through compareTwoAdults, it prints a confirmation message for the user 
        cout << "PASSED: Differentiated the better of two adults" << endl;
    }
    if (sortAdultVec(lotsOfCouts)){ //if it successfully makes it through sortAdultVec, it tells the user it was sucessful
      cout << "PASSED: Sorted a vector of adults as expected" << endl;
    } 
    else{
        cout << "Completed everything in sortAdultVec, but some of the tests FAILED" << endl;
    }
    if (runGeneticsUnitTests()){
        cout << "PASSED: All the genetics tests passed" << endl;
    }
    else{
        cout << "Unidentified error in testing_genetics" << endl;
    }
    // if(dominationTest()){
    //     cout << "dominationTest has ended" <<endl;
    // }
    //if(GDTest()){
     //   cout << "GDTest has ended" << endl;
    //}
    
    return 0;
}