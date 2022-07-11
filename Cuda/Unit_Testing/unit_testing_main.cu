#define UNITTEST
#include <iostream>
using std::cout;
using std::endl;
#include "../Unit_Testing/testing_sorts.h"
#include "../Unit_Testing/testing_genetics.h"
// nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.22.27905\bin\Hostx64\x64\cl.exe" -o unit_test unit_testing_main.cu -arch=compute_50 -code=sm_50 

int main(){
    //determines whether or much of the workings of the unit tests will be printed or not
    bool lotsOfCouts = true; 

    // Completes a series of comparisons on different pairs of Adults to ensure that rankSort and rankDistance sort are working as expected
    //      e.g. rankSort can identify which individual has a higher rank
    if (compareTwoAdults()){ 
        cout << "PASSED: Differentiated the better of two adults" << endl;
    }
    else{
        cout << "Could not correctly differentiate between two adults in all instances" << endl;
    }

    // Completes 6 test where a vector is sorted both by rankSort and rankDistanceSort
    // The results of these tests are compared to the expected order these vectors should appear in
    if (sortAdultVec(lotsOfCouts)){
      cout << "PASSED: Sorted a vector of adults as expected" << endl;
    } 
    else{
        cout << "Completed everything in sortAdultVec, but some of the tests FAILED" << endl;
    }

    // Calls a function that initializes a set of cConstants for unit testing and completes tests to verify
    //      different aspects of the genetic algorithm, ga_crossover, and annealing
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