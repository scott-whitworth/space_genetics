#include "../Genetic_Algorithm/adult.h"
//#include UNITTEST
#include <iostream>
#include <vector>
#include <string> //allows us to use to_string and string
using std::cout;
using std::endl;


// this function will check if giveRank is functioning as expected
bool dominationTest();

//test version of the giveRank function from optimization
void giveRankTest(std::vector<Adult> & allAdults);

// this function will run through giveDistanceTest to see if it is working correctly
bool GDTest();

//test version of the giveDistance function from optimization
void giveDistanceTest(std::vector<Adult> pool, int poolSize);


bool dominationTest(){

    //vector for the test
    std::vector<Adult> GRtest;

    int r = 0;//rank
    int d = 1;//distance
    int vectSize = 8;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        GRtest.push_back(Adult(r, d)); 
    }
    
    //giving each adult specific posDiffs and speedDiffs for expected results
    //should be rank 5; expected dominatedByCount = 4
    GRtest[0].posDiff = 25;
    GRtest[0].speedDiff = 25;
    //should be rank 4; expected dominatedByCount = 2
    GRtest[1].posDiff = 15;
    GRtest[1].speedDiff = 15;
    //should be rank 3; expected dominatedByCount = 1?
    GRtest[2].posDiff = 1;
    GRtest[2].speedDiff = 1;
    //should be rank 2; expected dominatedByCount = 1
    GRtest[3].posDiff = 1;
    GRtest[3].speedDiff = 1;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[4].posDiff = 1;
    GRtest[4].speedDiff = 1;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[5].posDiff = 1;
    GRtest[5].speedDiff = 1;
    //
    GRtest[6].posDiff = 25;
    GRtest[6].speedDiff = 25;
    //
    GRtest[7].posDiff = 25;
    GRtest[7].speedDiff = 25;
    
    giveRankTest(GRtest);
    cout << "Ranks:" << endl;
    for(int i = 0; i < GRtest.size(); i++){
        cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << endl;  
    }
    
    //giveDistanceTest(GRtest, GRtest.size());


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
    domination.resize(allAdults.size());
    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    //Note: fill the vector with 0s to make sure the count is accurate
    //TODO: unit test to make sure the whole vector is actually initially filled with 0's and not just the first index or the original vector size
    std::vector<int> dominatedByCount; 
    std::vector<int> secondaryCount;

    for(int i = 0; i < allAdults.size(); i++){
        dominatedByCount.push_back(0);        
        //cout << dominatedByCount[i] << ", ";
    }

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){
        cout << "i: " << i << " ";
        
        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){
            cout << " j: " << j << " ";

            //Check to see if i dominates j
            if (dominationCheckTest(allAdults[i], allAdults[j])){
                //Put the jth index in the set of individuals dominated by i
                cout << " i dominated j ";
                domination[i].push_back(j);
                //cout << " done " << endl;
                //TODO: will this add too many to j's domination count? When it i's current value reaches j's current value it will have already recorded the dominaton here, but it will be recorded again 
                //Add one to j's dominated by count
                //dominatedByCount[j]++; 
                cout << "count: " << dominatedByCount[j] << endl;
            }
            //Check to see if j dominates i
            else if (dominationCheckTest(allAdults[j], allAdults[i])) {
                //TODO: this may have the same redundancy that was mentioned above with things being added to this vector too many times
                //Put the ith index in the set of individuals dominated by j
                //domination[j].push_back(i);
                //cout << " j dominated i ";
                //Add one to i's dominated by count
                dominatedByCount[i]++; 
                
            }
            
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            
            cout << " front size: " << front.size() << endl;
            front.push_back(i);
        }

        
        cout << " i was dominated by j " << dominatedByCount[i] << " times ";
        cout << endl;
    }

    //Used to assign rank number
    int rankNum = 1;
    //vector to store individuals' indexes in next front
    std::vector<int> newFront;
    cout << " front size: " << front.size();
    cout << " domination size: " << domination.size() << endl;

    for(int k = 0; k < domination.size(); k++){
        for(int l = 0; l < domination[k].size(); l++){
            cout << domination[k][l] << ", ";
        } 
        cout << endl;
    }
    //int p = 0;

    for(int i = 0; i < dominatedByCount.size(); i++){
            secondaryCount.push_back(dominatedByCount[i]);
            cout << dominatedByCount[i] << ", ";
    }
    cout << endl;
    for(int i = 0; i < secondaryCount.size(); i++){
            cout << secondaryCount[i] << ", ";
    }
    cout << endl;
    //go until all individuals have been put in better ranks and none are left to give a ranking
    while(!front.empty()) {
        //empty the new front to put new individuals in
        std::vector<int>().swap(newFront);
        
        //loop through all individuals in old front
        //These individuals already have their rank set
        for(int k = 0; k < front.size(); k++){
            
            cout << " k: " << k;

            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < domination[k].size(); l++){
                cout << " l: " << l;

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[k][l]]--;
                cout << " index: " << domination[k][l] << " count: " << dominatedByCount[domination[k][l]] << " ";
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[k][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[k][l]].rank = rankNum + 1;
                    cout << " index: " << domination[k][l] << " current rank: " << allAdults[domination[k][l]].rank << " ";
                    //Add the index of the this adult to newFront

                    newFront.push_back(l);
                }
                
            }
            cout << endl;
        }

        //increment the rank number
        if(!newFront.empty()){
            rankNum++;
        }
        cout << " current overall rank: " << rankNum << endl;
        //empty the current front
        std::vector<int>().swap(front);

        std::vector<int>().swap(dominatedByCount);

        
        for(int i = 0; i < secondaryCount.size(); i++){
            if(secondaryCount[i] >= 0){
                secondaryCount[i]--;
            }
            dominatedByCount.push_back(secondaryCount[i]);
        }

        cout << "new counts: " << endl;
        for(int i = 0; i < dominatedByCount.size(); i++){
            cout << dominatedByCount[i] << ", ";
        }
        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        
        front = newFront;
        cout << endl;
        
        if(front.empty()){ 
            for(int i = 0; i < dominatedByCount.size(); i++){
                if(dominatedByCount[i] > 0){
                    
                    front.resize(front.size() + 1);
                }
            }    
        }
        
    }
}


bool GDTest(){

    //vector for the test
    std::vector<Adult> testGD;

    int r = 1;//rank
    int d = 1;//distance
    int vectSize = 6;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        testGD.push_back(Adult(r, d));
        testGD[i].posDiff = i + 5;
        testGD[i].speedDiff = i + 5;
    }
    
    giveDistanceTest(testGD, testGD.size());
    cout << "Distance:" << endl;
    for(int i = 0; i < testGD.size(); i++){
        cout << "d: " << testGD[i].distance << endl;  
    }

    return true;
}


//gives each adult in the pool a distance representing how different it is from other individuals
void giveDistanceTest(std::vector<Adult> pool, int poolSize){

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(pool.begin(), pool.begin() + poolSize, rankSort);

    for (int i = 0; i < poolSize; i++ ){
        //reset each individual's distance
        pool[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(pool.begin(), pool.begin() + poolSize, LowerPosDiff);
    //Set the boundaries
    pool[0].distance = 1e+12;
    pool[poolSize - 1].distance = 1e+12;


    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    for(int i = 1; i < poolSize - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffLeft = pool[i+1].posDiff/pool[poolSize - 1].posDiff;
        normalPosDiffRight = pool[i-1].posDiff/pool[poolSize - 1].posDiff;
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(pool[poolSize - 1].posDiff - pool[0].posDiff));
        cout << i << " posDiff Distance " << pool[i].distance << endl;
    }

    //Repeat above process for speedDiff    
    std::sort(pool.begin(), pool.begin() + poolSize, LowerSpeedDiff);
    //Set the boundaries
    pool[0].distance = 1e+12;
    pool[poolSize - 1].distance = 1e+12;

    
    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalSpeedDiffLeft;
    double normalSpeedDiffRight;
    for(int i = 1; i < poolSize - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalSpeedDiffLeft = pool[i+1].speedDiff/pool[poolSize - 1].speedDiff;
        normalSpeedDiffRight = pool[i-1].speedDiff/pool[poolSize - 1].speedDiff;
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(pool[poolSize - 1].speedDiff - pool[0].speedDiff));

        cout << i << " total Distance " << pool[i].distance << endl;
    }
}