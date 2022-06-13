#define UNITTEST
#include "../Unit_Testing/testing_sorts.h"
#include "../Unit_Testing/testing_genetics.h"

int main(){
    bool lotsOfCouts = false; //a boolean variable that tells determines whether or not all adults' ranks, distances, and statuses will be printed
    if (compareTwoAdults()){ //if it sucessfully makes it through compareTwoAdults, it prints a confirmation message for the user 
        cout << "PASSED: Differentiated the better of two adults" << endl;
    }
    if (sortAdultVec(lotsOfCouts)){ //if it successfully makes it through sortAdultVec, it tells the user it was sucessful
      cout << "PASSED: Sorted a vector of adults as expected" << endl;
    }
    if (firstParentsTest()){
        cout << "PASSED: Children can be converted to adults that can be sorted" << endl;
    }
    if(dominationTest()){
        cout << "dominationTest has ended" <<endl;
    }
    //if(GDTest()){
     //   cout << "GDTest has ended" << endl;
    //}
    
    return 0;
}