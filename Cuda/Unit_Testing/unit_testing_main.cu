#define UNITTEST
#include <iostream>
using std::cout;
using std::endl;
#include "../Unit_Testing/testing_sorts.h"
#include "../Unit_Testing/testing_genetics.h"

int main(){
    //determines whether or much of the workings of the unit tests will be printed or not
    bool lotsOfCouts = true; 
    if (compareTwoAdults()){ //if it sucessfully makes it through compareTwoAdults, it prints a confirmation message for the user 
    // TODO: The above comment, while correct, is not as helpful as something like:
    //       Check: <parameters of the unit test (I still don't really know exactly what this is testing)>
        cout << "PASSED: Differentiated the better of two adults" << endl;
    }
    else{
        cout << "Could not correctly differentiate between two adults in all instances" << endl;
    }
    if (sortAdultVec(lotsOfCouts)){ //if it successfully makes it through sortAdultVec, it tells the user it was sucessful
    //TODO: Same thing... what does 'successfully' mean here? We can see 'tell the user' based on the following cout 
      cout << "PASSED: Sorted a vector of adults as expected" << endl;
    } 
    else{
        cout << "Completed everything in sortAdultVec, but some of the tests FAILED" << endl;
    }
    if (runGeneticsUnitTests(lotsOfCouts)){
        cout << "PASSED: All the genetics tests passed" << endl;
    }
    else{
        cout << "Completed all the genetics units tests, but some FAILED" << endl;
    }
    // if(dominationTest()){
    //     cout << "dominationTest has ended" <<endl;
    // }
    //if(GDTest()){
     //   cout << "GDTest has ended" << endl;
    //}
    
    return 0;
}