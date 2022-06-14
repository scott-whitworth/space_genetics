#include "../Genetic_Algorithm/adult.h"
#include "../Genetic_Algorithm/individuals.h"
//#include UNITTEST
#include <iostream>
#include <iomanip>
#include <vector>
#include <string> //allows us to use to_string and string
#include <time.h>
using std::cout;
using std::endl;


// this function will check if giveRank is functioning as expected
//Input: None, but the class and posDiff and speedDiff of each adult will be created here
//Output: The completion of the test from calling the giveRankTest and giveDistanceTest
bool dominationTest();

//test version of the giveRank function from optimization, does not use cConstants
//Input: The adult vector created in dominationTest
//Output: Ranks given to the adults as well as cout statements
void giveRankTest(std::vector<Adult> & allAdults);

//test version of the giveDistance function from optimization, does not use cConstants
//Input: adult vetor from dominationTest and its size
//Output: distances given to each adult based on its speed and position compared to the adults around them
void giveDistanceTest(std::vector<Adult> & pool, int poolSize);

//giveRank fucntion that worked using arrays from individual, does not use cConstants
//Input: a set of individuals and its size
//Output: ranks for each individual
void oldGiveRankTest(Individual * pool, int poolSize);

//for testing the changeAnneal and changeInBest functions
//Input: None, but uses the same adult class setup as domination test and manually sets posDiff and speedDiff
//Output: Shows how the changeAnneal is working 
bool CAtest();

//test version of changeAnneal from optimization.cu, does not use cConstants
//Input: uses the adult vector and all the anneal information needed, 
//       generation is used for adjusting the anneal further every set amount of generations
//       previousBestPosDiff is found from the changeInBestTest and used to see if the best adult has changed or not
//       posTolerance is used for adjusting the annealing
//       dRate is the distinguishable rate which is adjusted based on the current best posDiff -- TODO: this part used to be cost based, what should it be now?
void changeAnnealTest(std::vector<Adult> newAdults, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate);

//test version of change in best, just checks if there has been a chnage in the best posDiff(used to be cost)
//Input: previousBestPosDiff (was cost) 
//       currentBest (the best posDiff of this generation)
//       dRate to compare to current and previous best
bool changeInBestTest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate);



bool dominationTest(){

    //vector for the test
    std::vector<Adult> GRtest;

    int r = 0;//rank
    double d = 1.0;//distance
    int test = 4;// 1: set of 6 different posDiffs and speedDiffs that leads to a normal expected outcome
                 // 2: set of 6 different posDiffs and speedDiffs with a pair of nans, this messes up giveDistanceTest
                 // 3: set of 1000 random posDiffs and speedDiffs of values 1-1000, meant to show giveRank works for large sets
                 // 4: set of 6 meant for testing the giveDistanceTest and seeing of it is giving the correct distances

    //creating the vector and filling it with adults
    //test 1 is used for showing a good example of how giveRank works well with kind values
    //test 1 is also meant for showing an area where giveDistance doesn't quite work well
    if(test == 4){
        int vectSize = 6;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }

        //Rank: 1
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 4;
        //Rank: 1
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;
        //Rank: 2
        GRtest[2].posDiff = 5;
        GRtest[2].speedDiff = 7;
        //Rank: 3
        GRtest[3].posDiff = 9.72;
        GRtest[3].speedDiff = 9.72;
        //Rank: 2
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;
        //Rank: 1
        GRtest[5].posDiff = 3;
        GRtest[5].speedDiff = 3;

    }else if (test == 2){//nans are used in this test to show how both giveRank and giveDistance fail
        int vectSize = 6;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }

        //Rank: 1
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 4;
        //Rank: 1
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;
        //Rank: 1?
        GRtest[2].posDiff = nan("");
        GRtest[2].speedDiff = nan("");
        //Rank: 3
        GRtest[3].posDiff = 9.72;
        GRtest[3].speedDiff = 9.72;
        //Rank: 2
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;
        //Rank: 1
        GRtest[5].posDiff = 3;
        GRtest[5].speedDiff = 3;
    }else if(test == 3){//larger test with random values
        int vectSize = 1000;//size of the vector for the test
        srand(time(NULL));//random seed
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
            //random values from 1 to 1000
            GRtest[i].posDiff = 1 + (rand() % 999);
            GRtest[i].speedDiff = 1 +  (rand() % 999);
        }
        //GRtest[100].posDiff = nan("");
        //GRtest[100].speedDiff = nan("");
        //GRtest[200].posDiff = nan("");
        //GRtest[200].speedDiff = nan("");
    }else if (test == 4){
        int vectSize = 6;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }
        //Rank: 1
        GRtest[0].posDiff = 1;
        GRtest[0].speedDiff = 4;
        //GRtest[0].distance = 3;
        //Rank: 1
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;
        //GRtest[1].distance = 3;
        //Rank: 2
        GRtest[2].posDiff = 5;
        GRtest[2].speedDiff = 7;
        //Rank: 3
        GRtest[3].posDiff = 9.72;
        GRtest[3].speedDiff = 9.72;
        //Rank: 2
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;
        //Rank: 1
        GRtest[5].posDiff = 3;
        GRtest[5].speedDiff = 3;
    }
    else{
        cout << "Not a valid test" << endl;
    }

    //call the giveRankTest
    giveRankTest(GRtest);

    //print the results
    std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankSort);
    cout << "Ranks:" << endl;
    for(int i = 0; i < GRtest.size(); i++){
        if(GRtest.size() > 50){//print every 50th adult for large vector
            if(i % 50 == 0){
                cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << endl;
            }
        }else{
           cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << endl; 
        }
    }

    
    //get the distance for each adult now that the ranks are set
    giveDistanceTest(GRtest, GRtest.size());
    std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankDistanceSort);
    //print the results
    cout << endl << "Total Distance: " << endl;
    for(int i = 0; i < GRtest.size(); i++){
        if(GRtest.size() > 50){//print every 50th adult for large vector
            if(i % 50 == 0){
                cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << " d:" << GRtest[i].distance << endl;
            }
        }else{
            cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << " d:" << GRtest[i].distance << endl;
        }
    }
    
    //same size as whatever test is running
    int vectSize = GRtest.size();
    //edited individual class
    Individual * testPool = new Individual[GRtest.size()];

    if(test == 3){//larger test with random values
        srand(time(NULL));
        for(int i = 0; i < vectSize; i++){
            //edited individual class
            testPool[i] = Individual(r, d);
            //random values from 1 to 1000, same as adult class
            testPool[i].posDiff = GRtest[i].posDiff;
            testPool[i].speedDiff = GRtest[i].speedDiff;
        }
    }else if(test == 2){
        for(int i = 0; i < vectSize; i++){
            //edited individual class
            testPool[i] = Individual(r, d);
            //same to compare ranks
            testPool[i].posDiff = GRtest[i].posDiff;
            testPool[i].speedDiff = GRtest[i].speedDiff;
        }
    }else if (test == 1){
        for(int i = 0; i < vectSize; i++){
            testPool[i] = Individual(r, d);
            //same to compare ranks
            testPool[i].posDiff = GRtest[i].posDiff;
            testPool[i].speedDiff = GRtest[i].speedDiff;
        }
    }else if (test == 4){
        for(int i = 0; i < vectSize; i++){
            testPool[i] = Individual(r, d);
            //same to compare ranks
            testPool[i].posDiff = GRtest[i].posDiff;
            testPool[i].speedDiff = GRtest[i].speedDiff;
        }
    }

    //do this sort after so the correct indexes are copied to testPool
    std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankSort);
    cout << "Best-> Rank: " << GRtest[0].rank << ", p: " << GRtest[0].posDiff << ", s: " << GRtest[0].speedDiff << ", d: " << GRtest[0].distance << endl;
    cout << "2nd Worst-> Rank: " << GRtest[GRtest.size() - 2].rank << ", p: " << GRtest[GRtest.size() - 2].posDiff << ", s: " << GRtest[GRtest.size() - 2].speedDiff << ", d: " << GRtest[GRtest.size() -2].distance << endl;
    cout << "Worst-> Rank: " << GRtest[GRtest.size() - 1].rank << ", p: " << GRtest[GRtest.size() - 1].posDiff << ", s: " << GRtest[GRtest.size() - 1].speedDiff << ", d: " << GRtest[GRtest.size() -1].distance << endl;

    //call the old giveRank function that uses arrays
    oldGiveRankTest(testPool, vectSize);
    //print the results
    cout << "Ranks:" << endl;
    for(int i = 0; i < vectSize; i++){
        if(vectSize > 50){
            if(i % 50 == 0){
                cout << "r: " << testPool[i].rank << ", p: " << testPool[i].posDiff << ", s: " << testPool[i].speedDiff << endl;
            }
        }else{
            cout << "r: " << testPool[i].rank << ", p: " << testPool[i].posDiff << ", s: " << testPool[i].speedDiff << endl;
        }
        
    }
    //compare the best and worst and how they were ranked
    std::sort(testPool, testPool + vectSize, rankSortInd);
    cout << "Best-> Rank: " << testPool[0].rank << ", p: " << testPool[0].posDiff << ", s: " << testPool[0].speedDiff << endl;
    cout << "2nd Worst-> Rank: " << testPool[GRtest.size() - 2].rank << ", p: " << testPool[GRtest.size() - 2].posDiff << ", s: " << testPool[GRtest.size() - 2].speedDiff << endl;
    cout << "Worst-> Rank: " << testPool[GRtest.size() - 1].rank << ", p: " << testPool[GRtest.size() - 1].posDiff << ", s: " << testPool[GRtest.size() - 1].speedDiff << endl;
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
    //crashed when the size was not set for both dimensions
    domination.resize(allAdults.size());
    
    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    //Note: fill the vector with 0s to make sure the count is accurate
    //TODO: unit test to make sure the whole vector is actually initially filled with 0's and not just the first index or the original vector size
    std::vector<int> dominatedByCount; 
    //crashed when this was not resized and filled with zeros (cant add 1 to an index that doesn't exist)
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){

            //Check to see if i dominates j
            if (dominationCheckTest(allAdults[i], allAdults[j]) || isnan(allAdults[j].posDiff) || isnan(allAdults[j].speedDiff)){
                //Put the jth index in the set of individuals dominated by i

                domination[i].push_back(j);
                
                

                //domination[i][j] = j;//this does not work
            }
            //Check to see if j dominates i
            else if (dominationCheckTest(allAdults[j], allAdults[i])  || isnan(allAdults[i].posDiff) || isnan(allAdults[i].speedDiff)){
                
                dominatedByCount[i]++;//tried it with this in the first if statement and it does not work, leads to a miscount
                
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            //if(isnan(allAdults[i].posDiff)  || isnan(allAdults[i].speedDiff)){
                //dominatedByCount[i] = allAdults.size()-1;
            //}else{
                allAdults[i].rank = 1;
                front.push_back(i);
            //}
            
        }
    }

    //Used to assign rank number
    int rankNum = 1;
    //vector to store individuals' indexes in next front
    std::vector<int> newFront;
    if(allAdults.size() > 20){
        cout << "***Done with domination process***" << endl;
        cout << " front size: " << front.size();
        cout << " domination size: " << domination.size() << endl;
    }else{
        cout << " front size: " << front.size();
        cout << " domination size: " << domination.size() << endl;
        //print the 2d vector that has the index of every adult that has been dominated at least once
        for(int k = 0; k < domination.size(); k++){
            for(int l = 0; l < domination[k].size(); l++){
                cout << domination[k][l] << ", ";
            } 
        cout << endl;
        }
        //print the dominatedByCount of every adult
        for(int i = 0; i < dominatedByCount.size(); i++){
            cout << dominatedByCount[i] << ", ";
        }
    }
    

    
    //go until all individuals have been put in better ranks and none are left to give a ranking
    while(!front.empty()) {
        //empty the new front to put new individuals in
        std::vector<int>().swap(newFront);
        
        //loop through all individuals in old front
        //These individuals already have their rank set
        for(int k = 0; k < front.size(); k++){
            
            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < domination[front[k]].size(); l++){
                
                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[front[k]][l]]--;
                
                //dominatedByCount[domination[front[k]][l]]--;
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[front[k]][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[front[k]][l]].rank = rankNum + 1;
                    //Add the index of the this adult to newFront
                    //TODO: these counts are still getting into the negatives
                    newFront.push_back(domination[front[k]][l]);
                }
                
            }
        }

        //increment the rank number
        rankNum++;
        //cout << " current overall rank: " << rankNum << endl;
        //empty the current front
        std::vector<int>().swap(front);

        //print the counts to check for negatives
        if(allAdults.size() > 20){
            //no cout
        }else{
            cout << "new counts: " << endl;
            for(int i = 0; i < dominatedByCount.size(); i++){
                cout << dominatedByCount[i] << ", ";
            }
        }
        
        //print the front to see what indexes are in it
        if(allAdults.size() > 20){
            //no cout
        }else{
            cout << endl << "new front: " << endl;
            for(int i = 0; i < newFront.size(); i++){
                cout << newFront[i] << ", ";
            }
        }
        
        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
        //cout << endl;
    }
    cout << endl << "Index 0 dominated count: " << dominatedByCount[0] << endl;
    cout << endl << "Index 1 dominated count: " << dominatedByCount[1] << endl;
}


//gives each adult in the pool a distance representing how different it is from other individuals
void giveDistanceTest(std::vector<Adult> & pool, int poolSize){

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
    //could we instead have it find the distance from i = 0 and i = 1
    for(int i = 1; i < poolSize-1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffRight = pool[i-1].posDiff/pool[poolSize - 1].posDiff;
        normalPosDiffLeft = pool[i+1].posDiff/pool[poolSize - 1].posDiff;
        
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(pool[poolSize - 1].posDiff - pool[0].posDiff));
        //cout << i << " posDiff Distance " <<  pool[i].distance << endl;
    }

    //Repeat above process for speedDiff    
    std::sort(pool.begin(), pool.begin() + poolSize, LowerSpeedDiff);
    //Set the boundaries
    pool[0].distance = 1e+12;
    pool[poolSize - 1].distance = 1e+12;
    //cout << endl << "Post lowerSpeedDiff sort:" << endl;
    for(int i = 0; i < pool.size(); i++){
        if(i % 1 == 0){
           //cout << "r: " << pool[i].rank << ", p: " << pool[i].posDiff << ", s: " << pool[i].speedDiff << " d:" << pool[i].distance << endl; 
        }
          
    }

    
    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalSpeedDiffLeft;
    double normalSpeedDiffRight;
    for(int i = 1; i < poolSize-1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalSpeedDiffRight = pool[i-1].speedDiff/pool[poolSize - 1].speedDiff;
        normalSpeedDiffLeft = pool[i+1].speedDiff/pool[poolSize - 1].speedDiff;
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(pool[poolSize - 1].speedDiff - pool[0].speedDiff));

    }

    
}

void oldGiveRankTest(Individual * pool, int poolSize) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf
    //Used to store the current front of individuals. first filled with the first front individuals(best out of all population)
    // filled with index of individuals in pool
    std::vector<int> front;
    
    //loop through each individual
    for (int i = 0; i < poolSize; i++){
        
        //number of times pool[i] has been dominated
        pool[i].dominatedCount = 0;
        //set of solutions that pool[i] dominates. Need to empty for each generation
        std::vector<int>().swap(pool[i].dominated);
        for(int j = 0; j < poolSize; j++){
            
            //if i dominates j, put the j index in the set of individuals dominated by i.
            if (dominationCheckTestv2(pool[i], pool[j])  || isnan(pool[j].posDiff) || isnan(pool[j].speedDiff)){
                pool[i].dominated.push_back(j);
            }
            //if j dominates i, increase the number of times that i has been dominated
            else if (dominationCheckTestv2(pool[j], pool[i])  || isnan(pool[i].posDiff) || isnan(pool[j].speedDiff)) {
                pool[i].dominatedCount++;
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (pool[i].dominatedCount == 0){
            pool[i].rank = 1;
            front.push_back(i);
        }
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
        for(int k = 0; k < front.size(); k++){
            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < pool[front[k]].dominated.size(); l++){
                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                pool[pool[front[k]].dominated[l]].dominatedCount--;
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (pool[pool[front[k]].dominated[l]].dominatedCount == 0){
                    pool[pool[front[k]].dominated[l]].rank = rankNum + 1;
                    newFront.push_back(pool[front[k]].dominated[l]);                        
                }
            }
        }
        //increment the rank number
        rankNum++;
        //empty the current front
        std::vector<int>().swap(front);
        //go to next front
        front = newFront;
    }
    cout << endl << "Index 0 dominated count (old): " << pool[0].dominatedCount << endl;
}

bool CAtest(){
    std::vector<Adult> CAtest;

    int r = 0;//rank
    double d = 1.0;//distance
    int vectSize = 6;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        CAtest.push_back(Adult(r, d)); 
        //GRtest[i].posDiff = i;
        //GRtest[i].speedDiff = i;
    }

    //Rank: 1
    CAtest[0].posDiff = 1e-6;
    CAtest[0].speedDiff = 1;
    //Rank: 1
    CAtest[1].posDiff = 4;
    CAtest[1].speedDiff = 1;
    //Rank: 2
    CAtest[2].posDiff = 5;
    CAtest[2].speedDiff = 7;
    //Rank: 3
    CAtest[3].posDiff = 9.72;
    CAtest[3].speedDiff = 9.72;
    //Rank: 2
    CAtest[4].posDiff = 11;
    CAtest[4].speedDiff = 2;
    //Rank: 1
    CAtest[5].posDiff = 3;
    CAtest[5].speedDiff = 3;
    //
    //GRtest[6].posDiff = 10;
    //GRtest[6].speedDiff = 10;
    //
    //GRtest[7].posDiff = 5;
    //GRtest[7].speedDiff = 5;
    //give them a rank
    giveRankTest(CAtest);

    //give them a distance (this distance will not be used?)
    giveDistanceTest(CAtest, CAtest.size());

    //**These are all based on what is in the config/cConstants**
    double anneal_initial = 0.10;
    double new_anneal;
    double currentAnneal = anneal_initial;
    double anneal_min = anneal_initial;
    double previousBestPosDiff = 0;
    double generation = 1;
    double posTolerance = 1e-14;
    double dRate = 1e-8;

    while(generation < 10){
        //call the test to adjust the anneal
        changeAnnealTest(CAtest, new_anneal, currentAnneal, anneal_min, previousBestPosDiff, generation, posTolerance, dRate);
        //print the results
        cout << "done with run " << generation << endl;
        cout << " new_anneal: " << new_anneal << " currentAnneal: " << currentAnneal << " anneal_min: " << anneal_min << endl;
        cout <<  " previousBestPosDiff: " << previousBestPosDiff << " dRate: " << dRate <<  endl;
        //adjust the best posDiff every generation to get more changes in the anneal
        CAtest[0].posDiff = CAtest[0].posDiff/10;
        //increase the generation 
        generation++;
    }

    return true;

}


void changeAnnealTest(std::vector<Adult> newAdults, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate){
    // Scaling anneal based on proximity to tolerance
    //cudaConstants (all based on what is in cConstants/config):
    int missionType = Rendezvous;
    double anneal_final = 1.0e-7;
    int change_check = 1;
    double anneal_factor = 0.75;
    double anneal_initial = 0.10;


    // Far away: larger anneal scale, close: smaller anneal
    if (missionType == Impact) {
        //Impact is only based on posDiff, so proximity-based annealing only relies on how close posDiff is to tolerance.
        new_anneal = currentAnneal * (1 - (posTolerance / newAdults[0].posDiff));
    }

    else if (missionType == Rendezvous) {
        if (posTolerance < newAdults[0].posDiff){ 
            //TO DO: decide what we want to do with this annealing   
            //Exponentially changing annealing, as oppose to what?
            // this may be where the anneal is jumping?
            new_anneal = currentAnneal * (1 - pow(posTolerance / newAdults[0].posDiff,2.0));
            if (new_anneal < anneal_final){
                new_anneal = anneal_final; //Set a true minimum for annealing
            }
        }
    }
    
    //Process to see if anneal needs to be adjusted
    // If generations are stale, anneal drops
    Adult currentBest;
    // Compare current best individual to that from CHANGE_CHECK (50) many generations ago.
    // If they are the same, change size of mutations
    if (static_cast<int>(generation) % (change_check) == 0) { 
        currentBest = newAdults[0];
        // checks for anneal to change
        // previousBest starts at 0 to ensure changeInBest = true on generation 0
        if ( !(changeInBestTest(previousBestPosDiff, currentBest, dRate)) ) { 
            //this ensures that changeInBest never compares two zeros, thus keeping dRate in relevance as the posDiff lowers
            if (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                while (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                    dRate = dRate/10; 
                }
                std::cout << "\nnew dRate: " << dRate << std::endl;
            }
            // If no change in BestIndividual across generations, reduce currentAnneal by anneal_factor while staying above anneal_min
            //reduce anneal_min
            // the exponential decrease here leads to a huge jump on the last generation?
            anneal_min = anneal_initial*exp(-sqrt(posTolerance/newAdults[0].posDiff)*generation);
            if (anneal_min < anneal_final){
                anneal_min = anneal_final;//Set a true minimum for annealing
            }

            //Rendezvous mission uses anneal_min, Impact does not
            if(missionType == Impact) {
                currentAnneal = currentAnneal * anneal_factor;
            }
            else if (missionType == Rendezvous){
                currentAnneal = (currentAnneal * anneal_factor > anneal_min)? (currentAnneal * anneal_factor):(anneal_min);
            }
            std::cout << "\nnew anneal: " << currentAnneal << std::endl;              
        }

        previousBestPosDiff = currentBest.posDiff; //posDiff here used to be cost
    }
}

bool changeInBestTest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate) {
    //truncate is used here to compare doubles via the distinguguishRate, to ensure that there has been relatively no change.
        if (trunc(previousBestPosDiff/distinguishRate) != trunc(currentBest.posDiff/distinguishRate)) {
            cout << "true"<< endl;
            return true;
        }
        else { 
            cout << "false"<< endl;
            return false;
        }
}