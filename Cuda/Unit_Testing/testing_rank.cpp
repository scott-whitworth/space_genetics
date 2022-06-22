#include "../Genetic_Algorithm/adult.h"
//#include "../Genetic_Algorithm/individuals.h"
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
void giveRankTest(std::vector<Adult> & allAdults, int missionType);

//test version of the giveDistance function from optimization, does not use cConstants
//Input: adult vetor from dominationTest and its size
//Output: distances given to each adult based on its speed and position compared to the adults around them
void giveDistanceTest(std::vector<Adult> & pool, int poolSize, int missionType);

//giveRank fucntion that worked using arrays from individual, does not use cConstants
//Input: a set of individuals and its size
//Output: ranks for each individual
//void oldGiveRankTest(Individual * pool, int poolSize);

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
//NOTE: this header is outdated, this is the new changeAnneal function that uses cost
void changeAnnealTest(const std::vector<Adult>& oldAdults, int missionType, double & new_anneal, double & currentAnneal);

//function that calculates the cost
double calculateCostTest(const std::vector<Adult> & oldAdults, int missionType);

//test version of change in best, just checks if there has been a chnage in the best posDiff(used to be cost)
//Input: previousBestPosDiff (was cost) 
//       currentBest (the best posDiff of this generation)
//       dRate to compare to current and previous best
//bool changeInBestTest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate);
bool newDuplicateTest();


bool dominationTest(){

    //vector for the test
    std::vector<Adult> GRtest;

    int r = 0;//rank
    double d = 1.0;//distance
    int test = 6;// 1: set of 6 different posDiffs and speedDiffs that leads to a normal expected outcome
                 // 2: set of 6 different posDiffs and speedDiffs for the new version of dominationCheckTest that accounts for mission type
                 // 3: set of 1000 random posDiffs and speedDiffs of values 1-1000, meant to show giveRank works for large sets
                 // 4: error status test of 6 adults
                 // 5: set of easy values to calculate distance and test what happens when we sort differently
                 // 6: random version of test 6, goes through all 3 ways to calculate distance based on missionType
                 // 7: copied version of meg's problem set

    int missionType; //ren 1, imp 2

    //creating the vector and filling it with adults
    //test 1 is used for showing a good example of how giveRank works well with kind values
    //test 1 is also meant for showing an area where giveDistance doesn't quite work well
    if(test == 1){
        int vectSize = 6;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }

        //Rank (ren): 1
        //Rank (imp): 1
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 4;

        //Rank (ren): 1
        //Rank (imp): 3
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;

        //Rank (ren): 2
        //Rank (imp): 1
        GRtest[2].posDiff = 5;
        GRtest[2].speedDiff = 7;

        //Rank (ren): 3
        //Rank (imp): 1
        GRtest[3].posDiff = 9.72;
        GRtest[3].speedDiff = 9.72;

        //Rank (ren): 2
        //Rank (imp): 4
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;

        //Rank (ren): 1
        //Rank (imp): 2
        GRtest[5].posDiff = 3;
        GRtest[5].speedDiff = 3;

    }
    else if (test == 2){//this test is for the new version of dominationCheckTest that accounts for mission type
                        //this test used to be different from test 1, now it is the same
        int vectSize = 6;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }
        //higher speedDiffs will be better
        //Rank (ren): 1
        //Rank (imp): 1
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 4;

        //Rank (ren): 1
        //Rank (imp): 3
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;

        //Rank (ren): 2
        //Rank (imp): 1
        GRtest[2].posDiff = 5;
        GRtest[2].speedDiff = 7;

        //Rank (ren): 3
        //Rank (imp): 1
        GRtest[3].posDiff = 9.72;
        GRtest[3].speedDiff = 9.72;

        //Rank (ren): 2
        //Rank (imp): 4
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;

        //Rank (ren): 1
        //Rank (imp): 2
        GRtest[5].posDiff = 3;
        GRtest[5].speedDiff = 3;

    }
    else if(test == 3){//larger test with random values
        int vectSize = 1000;//size of the vector for the test
        srand(time(NULL));//random seed
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d, VALID)); 
            //random values from 1 to 1000
            GRtest[i].posDiff = 1 + (rand() % 1000);
            GRtest[i].speedDiff = 1 +  (rand() % 1000);
        }
        GRtest[100].errorStatus = SUN_ERROR;
        GRtest[200].errorStatus = SUN_ERROR;
        GRtest[300].errorStatus = SUN_ERROR;
        GRtest[400].errorStatus = SUN_ERROR;
    }
    else if(test == 4){// 4: error status test of 6 adults
        int vectSize = 8;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d, VALID)); 
        }

        //Rank: 4
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 4;
        GRtest[0].errorStatus = SUN_ERROR;
        //Rank: 4
        GRtest[1].posDiff = 4;
        GRtest[1].speedDiff = 1;
        GRtest[1].errorStatus = NAN_ERROR;
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
        //Rank: 1
        GRtest[6].posDiff = 2;
        GRtest[6].speedDiff = 4;
        //Rank: 5
        GRtest[7].posDiff = 10;
        GRtest[7].speedDiff = 12;
        GRtest[7].errorStatus = NAN_ERROR;
    }
    else if(test == 5){// 5: set of easy values to calculate distance and test what happens when we sort differently
        int vectSize = 8;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d)); 
        }

        //Rank (ren): 1
        //rank (imp): 3
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[0].posDiff = 2;
        GRtest[0].speedDiff = 1;

        //Rank (ren): 1
        //Rank (imp): 4
        //distance (imp low): 0.2625
        //distance (imp high): 1.2
        //distance (rendezvous): 0.2625
        GRtest[1].posDiff = 8;
        GRtest[1].speedDiff = 3;

        //Rank (ren): 2
        //Rank (imp): 1
        //distance (imp low): 0.775
        //distance (imp high): 10.15
        //distance (rendezvous): 0.775
        GRtest[2].posDiff = 9;
        GRtest[2].speedDiff = 12;

        //Rank (ren): 3
        //Rank (imp): 5
        //distance (imp low): 0.6375
        //distance (imp high): 3.45
        //distance (rendezvous): 0.6375
        GRtest[3].posDiff = 14;
        GRtest[3].speedDiff = 4;

        //Rank (ren): 2
        //Rank (imp): 6
        //distance (imp low): 0.2675
        //distance (imp high): 2.15
        //distance (rendezvous): 0.225
        GRtest[4].posDiff = 11;
        GRtest[4].speedDiff = 2;

        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 0.8
        //distance (imp high): 8.3
        //distance (rendezvous): 0.8
        GRtest[5].posDiff = 5;
        GRtest[5].speedDiff = 6;

        //Rank (ren): 2
        //Rank (imp): 5
        //distance (imp low): 0.1625
        //distance (imp high): 1.1
        //distance (rendezvous): 0.2125
        GRtest[6].posDiff = 11;
        GRtest[6].speedDiff = 3;

        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[7].posDiff = 20;
        GRtest[7].speedDiff = 16;

    }
    else if(test == 6){// 6: randomized and longer version of test 6, goes through all 3 ways to calculate distance based on missionType
        int vectSize = 30;//size of the vector for the test
        srand(time(NULL));//random seed
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(r, d, VALID)); 
            //random values from 1 to 100
            GRtest[i].posDiff = 1 + (rand() % 100);
            GRtest[i].speedDiff = 1 +  (rand() % 100);
        }

    }else if(test == 7){
        int vectSize = 10;//size of the vector for the test
        //fill a vector with the adults
        for(int i = 0; i < vectSize; i++){
            GRtest.push_back(Adult(200, d)); 
        }
        //Rank (ren): 1
        //rank (imp): 3
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[0].posDiff = 150;
        GRtest[0].speedDiff = 20;

        //Rank (ren): 1
        //Rank (imp): 4
        //distance (imp low): 0.2625
        //distance (imp high): 1.2
        //distance (rendezvous): 0.2625
        GRtest[1].posDiff = 120;
        GRtest[1].speedDiff = 90;

        //Rank (ren): 2
        //Rank (imp): 1
        //distance (imp low): 0.775
        //distance (imp high): 10.15
        //distance (rendezvous): 0.775
        GRtest[2].posDiff = 180;
        GRtest[2].speedDiff = 30;

        //Rank (ren): 3
        //Rank (imp): 5
        //distance (imp low): 0.6375
        //distance (imp high): 3.45
        //distance (rendezvous): 0.6375
        GRtest[3].posDiff = 20;
        GRtest[3].speedDiff = 70;

        //Rank (ren): 2
        //Rank (imp): 6
        //distance (imp low): 0.2675
        //distance (imp high): 2.15
        //distance (rendezvous): 0.225
        GRtest[4].posDiff = 110;
        GRtest[4].speedDiff = 40;

        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 0.8
        //distance (imp high): 8.3
        //distance (rendezvous): 0.8
        GRtest[5].posDiff = 30;
        GRtest[5].speedDiff = 90;

        //Rank (ren): 2
        //Rank (imp): 5
        //distance (imp low): 0.1625
        //distance (imp high): 1.1
        //distance (rendezvous): 0.2125
        GRtest[6].posDiff = 50;
        GRtest[6].speedDiff = 120;

        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[7].posDiff = 220;
        GRtest[7].speedDiff = 970;
        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[8].posDiff = 20;
        GRtest[8].speedDiff = 120;
        //Rank (ren): 1
        //Rank (imp): 2
        //distance (imp low): 1e+12
        //distance (imp high): 1e+12
        //distance (rendezvous): 1e+12
        GRtest[9].posDiff = 340;
        GRtest[9].speedDiff = 90;

    }
    else{
        cout << "Not a valid test" << endl;
    }

    
    //test 7 has special printing and calling needs
    if(test == 6 || test == 7){
        for(int run = 0; run < 3; run++){//run for the 3 cases
            if(run == 0){//normal rendezvous mission
                missionType = 1; //ren 1, imp 2
                //call giveRank to rank each adult
                giveRankTest(GRtest, missionType);

                //call giveDistance to get the distance of each adult
                giveDistanceTest(GRtest, GRtest.size(), missionType);

                //do a rankDistanceSort before printing
                std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankDistanceSort);
                //print the results
                cout << "Rendezvous mission: " << endl;
                for(int i = 0; i < GRtest.size(); i++){
                    cout << "Rank: " << GRtest[i].rank << " Distance: " << GRtest[i].distance << " PosDiff: " << GRtest[i].posDiff << " speedDiff: " << GRtest[i].speedDiff << " Index: " << i <<  endl;
                }
            }
            else if(run == 1){//impact mission with no speedDiff distance calculation
                missionType = 2; //ren 1, imp 2
                //call giveRank to rank each adult
                giveRankTest(GRtest, missionType);

                //call giveDistance to get the distance of each adult based only on posDiff
                giveDistanceTest(GRtest, GRtest.size(), missionType);
                //do a rankDistanceSort before printing

                std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankDistanceSort);
                cout << "Impact mission (no speedDiff distance calc): " << endl;
                for(int i = 0; i < GRtest.size(); i++){
                    cout << "Rank: " << GRtest[i].rank << " Distance: " << GRtest[i].distance << " PosDiff: " << GRtest[i].posDiff << " speedDiff: " << GRtest[i].speedDiff << " Index: " << i << endl;
                }

            }
            else{//impact mission with speedDiff calc
            /*
                missionType = 2; //ren 1, imp 2
                //call giveRank to rank each adult
                giveRankTest(GRtest, missionType);

                missionType = 1;//to find the distance based on posdiff and speedDiff using lowerPosDiff
                giveDistanceTest(GRtest, GRtest.size(), missionType);
                //do a rankDistanceSort before printing

                std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankDistanceSort);
                cout << "Impact mission (speedDiff calc): " << endl;
                for(int i = 0; i < GRtest.size(); i++){
                    cout << "Rank: " << GRtest[i].rank << " Distance: " << GRtest[i].distance << " PosDiff: " << GRtest[i].posDiff << " speedDiff: " << GRtest[i].speedDiff << " Index: " << i <<  endl;
                }
                */
            }
        }
    }
    else{//this will print for tests 1-6
        missionType = 2;//1 ren, 2 imp
        //call giveRank to rank each adult
        giveRankTest(GRtest, missionType);
        //only sort by rank
        std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankSort);
        cout << "Ranks:" << endl;
        for(int i = 0; i < GRtest.size(); i++){
            if(GRtest.size() > 50){//print every 50th adult for large vector
                if(i % 50 == 0){
                    cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << ", e: " << GRtest[i].errorStatus << endl;
                }
            }
            else{
                cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << ", e: " << GRtest[i].errorStatus << endl; 
            }
        }

        //get the distance for each adult now that the ranks are set
        giveDistanceTest(GRtest, GRtest.size(), missionType);
        std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankDistanceSort);
        //print the results
        cout << endl << "Total Distance: " << endl;
        for(int i = 0; i < GRtest.size(); i++){
            if(GRtest.size() > 50){//print every 50th adult for large vector
                if(i % 50 == 0){
                    cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << " d:" << GRtest[i].distance << ", e: " << GRtest[i].errorStatus << endl;
                }
            }
            else{
                cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << " d:" << GRtest[i].distance << ", e: " << GRtest[i].errorStatus << endl;
            }
        }
        //these are to verify that a nan/invalid adult goes to the worst rank
        std::sort(GRtest.begin(), GRtest.begin() + GRtest.size(), rankSort);
        cout << "Best-> Rank: " << GRtest[0].rank << ", p: " << GRtest[0].posDiff << ", s: " << GRtest[0].speedDiff << ", d: " << GRtest[0].distance << endl;
        cout << "2nd Worst-> Rank: " << GRtest[GRtest.size() - 2].rank << ", p: " << GRtest[GRtest.size() - 2].posDiff << ", s: " << GRtest[GRtest.size() - 2].speedDiff << ", d: " << GRtest[GRtest.size() -2].distance << endl;
        cout << "Worst-> Rank: " << GRtest[GRtest.size() - 1].rank << ", p: " << GRtest[GRtest.size() - 1].posDiff << ", s: " << GRtest[GRtest.size() - 1].speedDiff << ", d: " << GRtest[GRtest.size() -1].distance << endl;
    }
        
    return true;
    
}

void giveRankTest(std::vector<Adult> & allAdults, int missionType) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf
    
    //Used to store the current front of adults. first filled with the first front adults(best out of all population)
    // filled with index of adults in allAdults
    std::vector<int> front;


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
    std::vector<int> dominatedByCount; 
    //crashed when this was not resized and filled with zeros (cant add 1 to an index that doesn't exist)
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){
            
            if(allAdults[i].errorStatus != VALID && allAdults[j].errorStatus == VALID){
                dominatedByCount[i]++;//tried it with this in the first if statement and it does not work, leads to a miscount

            }
            else if(allAdults[j].errorStatus != VALID && allAdults[i].errorStatus == VALID){
                domination[i].push_back(j);
                
            }else if (dominationCheckTest(allAdults[i], allAdults[j], missionType)){//Check to see if i dominates j
                //Put the jth index in the set of individuals dominated by i

                domination[i].push_back(j);

                //domination[i][j] = j;//this does not work
            }
            //Check to see if j dominates i
            else if (dominationCheckTest(allAdults[j], allAdults[i], missionType)){
                
                dominatedByCount[i]++;//tried it with this in the first if statement and it does not work, leads to a miscount
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
                allAdults[i].rank = 1;
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
        //These individuals already have their rank set
        for(int k = 0; k < front.size(); k++){
            
            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < domination[front[k]].size(); l++){
                
                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[front[k]][l]]--;
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[front[k]][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[front[k]][l]].rank = rankNum + 1;
                    //Add the index of the this adult to newFront
                    //TODO: these counts are still getting into the negatives for smaller tests
                    newFront.push_back(domination[front[k]][l]);
                }
                
            }
        }

        //increment the rank number
        rankNum++;
        //cout << " current overall rank: " << rankNum << endl;
        //empty the current front
        std::vector<int>().swap(front);
        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
        //print the counts to check for negatives
        /*
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
        }*/
        
        
        
        //cout << endl;
    }
    //cout << endl << "Index 0 dominated count: " << dominatedByCount[0] << endl;
    //cout << endl << "Index 1 dominated count: " << dominatedByCount[1] << endl;
}


//gives each adult in the pool a distance representing how different it is from other individuals
void giveDistanceTest(std::vector<Adult> & allAdults, int poolSize, int missionType){

    std::vector<int> validPool;

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(allAdults.begin(), allAdults.begin() + poolSize, rankSort);

    for(int i = 0; i < allAdults.size(); i++){
        if(allAdults[i].errorStatus == VALID){
            validPool.push_back(i);
        }
    }
    

    for (int i = 0; i < poolSize; i++ ){
        //reset each individual's distance
        allAdults[i].distance = 0.0;
    }
    int newPoolSize;
    newPoolSize = validPool.size();
    //Sort by the first objective function, posDiff
    std::sort(allAdults.begin(), allAdults.begin() + newPoolSize, LowerPosDiff);
    //Set the boundaries
    allAdults[0].distance += 1;
    allAdults[newPoolSize - 1].distance += 1;
    
    //cout << endl << "Post lowerPosDiff sort:" << endl;
    for(int i = 0; i < allAdults.size(); i++){
        if(i % 1 == 0){
           //cout << "r: " << allAdults[i].rank << ", p: " << allAdults[i].posDiff << ", s: " << allAdults[i].speedDiff << " d:" << allAdults[i].distance << ", e: " << allAdults[i].errorStatus << endl; 
        }  
    }
    

    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    //could we instead have it find the distance from i = 0 and i = 1
    for(int i = 1; i < newPoolSize-1; i++){
        //Divide left and right individuals by the worst individual to normalize
        if(missionType == Impact){
            normalPosDiffRight = allAdults[i-1].posDiff/allAdults[newPoolSize - 1].posDiff;
            
            normalPosDiffLeft = allAdults[i+1].posDiff/allAdults[newPoolSize - 1].posDiff;
            //cout << "Right: " << normalPosDiffRight << " Left: " << normalPosDiffLeft << endl;
        }else{
            normalPosDiffRight = allAdults[i-1].posDiff/allAdults[newPoolSize - 1].posDiff;
            normalPosDiffLeft = allAdults[i+1].posDiff/allAdults[newPoolSize - 1].posDiff;
        }
        
        
        //distance = distance + abs((i+1) - (i-1))
        allAdults[i].distance = allAdults[i].distance + abs(normalPosDiffLeft - normalPosDiffRight);// /(pool[poolSize - 1].posDiff - pool[0].posDiff));
        //cout << i << " posDiff Distance " <<  pool[i].distance << endl;
    }

    //Repeat above process for speedDiff
    if(missionType == Impact){
        //cout << endl << "Post HigherSpeedDiff sort:" << endl;
        //std::sort(allAdults.begin(), allAdults.begin() + newPoolSize, HigherSpeedDiff);
    }else{
        //cout << endl << "Post lowerSpeedDiff sort:" << endl;
        std::sort(allAdults.begin(), allAdults.begin() + newPoolSize, LowerSpeedDiff);
        allAdults[0].distance += 1;
        allAdults[newPoolSize - 1].distance += 1;
    }
    
    //Set the boundaries
    
    
    for(int i = 0; i < allAdults.size(); i++){
        if(i % 1 == 0){
           //cout << "r: " << allAdults[i].rank << ", p: " << allAdults[i].posDiff << ", s: " << allAdults[i].speedDiff << " d:" << allAdults[i].distance << ", e: " << allAdults[i].errorStatus << endl; 
        }
          
    }
    
    
    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalSpeedDiffLeft;
    double normalSpeedDiffRight;
    for(int i = 1; i < newPoolSize-1; i++){
        //Divide left and right individuals by the worst individual to normalize
        if(missionType == Impact){
            //normalSpeedDiffRight = pool[newPoolSize - 1].speedDiff/pool[i-1].speedDiff;
            //normalSpeedDiffLeft = pool[newPoolSize - 1].speedDiff/pool[i+1].speedDiff;
        }else{
            normalSpeedDiffRight = allAdults[i-1].speedDiff/allAdults[newPoolSize - 1].speedDiff;
            normalSpeedDiffLeft = allAdults[i+1].speedDiff/allAdults[newPoolSize - 1].speedDiff;
            allAdults[i].distance = allAdults[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(pool[poolSize - 1].speedDiff - pool[0].speedDiff));
        }
        //distance = distance + abs((i+1) - (i-1))
        

    }

    //cout << "Pre sort: " << endl;
    for(int i = 0; i < allAdults.size(); i++){
        if(i % 1 == 0){
           //cout << "r: " << allAdults[i].rank << ", p: " << allAdults[i].posDiff << ", s: " << allAdults[i].speedDiff << " d:" << allAdults[i].distance << ", e: " << allAdults[i].errorStatus << endl; 
        }
    }
    
}

bool CAtest(){
    std::vector<Adult> CAtest;

    int r = 0;//rank
    double d = 1.0;//distance
    int vectSize = 1;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        CAtest.push_back(Adult(r, d)); 
        //GRtest[i].posDiff = i;
        //GRtest[i].speedDiff = i;
    }

    //Rank: 1
    CAtest[0].posDiff = 3.14159e-8;
    CAtest[0].speedDiff = 2.71828e-9;

    //give them a distance (this distance will not be used?)
    //giveDistanceTest(CAtest, CAtest.size());

    //**These are all based on what is in the config/cConstants**
    int missionType = 2;//1 for imp, 2 for ren
    double anneal_initial = 0.10;
    double new_anneal;
    double currentAnneal = anneal_initial;
    //double anneal_min = anneal_initial;
    //double previousBestPosDiff = 0;
    int generation = 1;
    //double posTolerance = 1e-14;
    //double dRate = 1e-8;

    while(calculateCostTest(CAtest, missionType) >= Impact){
        //call the test to adjust the anneal
        changeAnnealTest(CAtest, missionType, new_anneal, currentAnneal);
        //print the results
        if(generation %1==0){
        cout << "done with run " << generation << endl;
        cout << "currentAnneal: " << currentAnneal << endl;
        cout << "Current cost: " << calculateCostTest(CAtest, missionType) << endl;
        cout << "Best posDiff: " << CAtest[0].posDiff << ", Best speedDiff: " << CAtest[0].speedDiff << endl;
        }
        //adjust the best posDiff every generation to get more changes in the anneal
        CAtest[0].posDiff = CAtest[0].posDiff/1.2;
        CAtest[0].speedDiff = CAtest[0].speedDiff/1.2;
        //increase the generation 
        generation++;
    }
        changeAnnealTest(CAtest, missionType, new_anneal, currentAnneal);
        cout << "Total generations = " << generation << endl;
        cout << "currentAnneal: " << currentAnneal << endl;
        cout << "Current cost: " << calculateCostTest(CAtest, missionType) << endl;
        cout << "Best posDiff: " << CAtest[0].posDiff << ", Best speedDiff: " << CAtest[0].speedDiff << endl;

    return true;

}
void changeAnnealTest(const std::vector<Adult>& oldAdults, int missionType, double & new_anneal, double & currentAnneal){
    
    //Calculate the current cost for this generation
    double curCost = calculateCostTest(oldAdults, missionType);
    double anneal_initial = 0.1;
    double anneal_final = 1.0e-7;
    //Caluclate the new current anneal
    //It will be a linear decrease from the initial/max anneal based on how well the current cost is compared to the cost threshold
    //If we happen to be close to covergence, and the cost we get is -0.0001:
  //currentAnneal =               0.1          * (              1.000001                  )
    //if the signs were correct, this would still be only a small decrease when very close to convergence
//    currentAnneal = anneal_initial * (1 - std::abs(costThreshold / curCost));
//    currentAnneal = anneal_initial * (1 - pow(std::abs(costThreshold / curCost),2.0));
//    currentAnneal = anneal_initial * (1 - (1 / curCost));
    if(missionType == Rendezvous){
        currentAnneal = anneal_initial * (1 - pow((Rendezvous / curCost),4.0));
    }else if (missionType == Impact){
        currentAnneal = anneal_initial * (1 - (Impact / curCost));
    }
    
//    currentAnneal = anneal_initial * (1 - pow((Rendezvous / curCost),2.0));
    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < anneal_final)
    {
        currentAnneal = anneal_final;
    }
}

double calculateCostTest(const std::vector<Adult> & oldAdults, int missionType){
    //Create the cost double that will be returned
    double cost; 
    double pos_threshold = 1.0e-10;
    double speed_threshold = 1.0e-11;

    //How to calculate the cost will depend on the mission type
    //In general, cost will be (for each mission parameter, the difference/the goal) - the number of mission parameters
    //This will mean a cost of 0 or below signifies that the individual has hit the mission goals

    //For impacts, the cost only depends on posDiff
    //  NOTE: While the RD sorting favors high speed for impacts, convergence ultimately comes down to posDiff
    if (missionType == Impact) {
        //Calculate the cost based only on position
        //The goal is position threshold and the current status is the best individual's posDiff
    //  cost =      bigger number    /         1e-10              -   1 = large number - 1
    //  The run before convergence:
    //  cost =   number <= 1e-10     /         1e-10              -   1 = smallish number - 1
//        cost = (oldAdults[0].posDiff /pos_threshold) - Impact;
        cost = (oldAdults[0].posDiff /pos_threshold);
        if(cost < 1){
            cost = 1;
        }
        //this might better if we did:

    }
    //For rendezvous, the cost depends on both posDiff and speedDiff
    else {
        //Similarly to impact, calculate how far the best adult's speed and position diffs are away from the goal and subtract by the number of mission goals
//        cost = ((oldAdults[0].posDiff / pos_threshold) + (oldAdults[0].speedDiff / speed_threshold)) - Rendezvous;
        cost = ((oldAdults[0].posDiff / pos_threshold) + (oldAdults[0].speedDiff / speed_threshold));
        if(cost < 2){
            cost = 2;
        }
    }
    
    //Return the calculated cost
    return cost; 
}

bool newDuplicateTest(){
    //vector for the test
    std::vector<Adult> dupeTest;

    int vectSize = 8;//size of the vector for the test
    srand(time(NULL));//random seed
    //fill a vector with the adults
    for(int i = 0; i < vectSize; i++){
        //default constuctor, don't care about rank or distacne
        dupeTest.push_back(Adult()); 
        //random values from 1 to 100
        //dupeTest[i].posDiff = 1 + (rand() % 100);
        //dupeTest[i].speedDiff = 1 +  (rand() % 100);
    }
    //Rank (ren): 1
    //rank (imp): 3
    //distance (imp low): 1e+12
    //distance (imp high): 1e+12
    //distance (rendezvous): 1e+12
    dupeTest[0].posDiff = 5;
    dupeTest[0].speedDiff = 6;

    //Rank (ren): 1
    //Rank (imp): 4
    //distance (imp low): 0.2625
    //distance (imp high): 1.2
    //distance (rendezvous): 0.2625
    dupeTest[1].posDiff = 5;
    dupeTest[1].speedDiff = 6;

    //Rank (ren): 2
    //Rank (imp): 1
    //distance (imp low): 0.775
    //distance (imp high): 10.15
    //distance (rendezvous): 0.775
    dupeTest[2].posDiff = 5;
    dupeTest[2].speedDiff = 6;

    //Rank (ren): 3
    //Rank (imp): 5
    //distance (imp low): 0.6375
    //distance (imp high): 3.45
    //distance (rendezvous): 0.6375
    dupeTest[3].posDiff = 5;
    dupeTest[3].speedDiff = 6;

    //Rank (ren): 2
    //Rank (imp): 6
    //distance (imp low): 0.2675
    //distance (imp high): 2.15
    //distance (rendezvous): 0.225
    dupeTest[4].posDiff = 11;
    dupeTest[4].speedDiff = 2;

    //Rank (ren): 1
    //Rank (imp): 2
    //distance (imp low): 0.8
    //distance (imp high): 8.3
    //distance (rendezvous): 0.8
    dupeTest[5].posDiff = 5;
    dupeTest[5].speedDiff = 6;

    //Rank (ren): 2
    //Rank (imp): 5
    //distance (imp low): 0.1625
    //distance (imp high): 1.1
    //distance (rendezvous): 0.2125
    dupeTest[6].posDiff = 11;
    dupeTest[6].speedDiff = 2;

    //Rank (ren): 1
    //Rank (imp): 2
    //distance (imp low): 1e+12
    //distance (imp high): 1e+12
    //distance (rendezvous): 1e+12
    dupeTest[7].posDiff = 5;
    dupeTest[7].speedDiff = 6;

    //Vector that will hold the adults who are potential parents
    //The criteria for being a parent is being in the top survivor_count number of adults in the oldAdult pool and not being a duplicate (distance > 0)
    //      Duplicate adults will generate children in a separate fashion
    std::vector<Adult> parents; 

    //Vector for duplicates, based on the same criteria as above
    std::vector<Adult> duplicates;
    //this sort is not nescessary
    std::sort(dupeTest.begin(), dupeTest.end(), rankDistanceSort);
    for (int i = 0; i < vectSize; i++){//loop through all the adults
        for(int j = i+1; j < vectSize; j++){//i+1 so it doesn't check itself or past indexes
            if(duplicateCheckTest(dupeTest[i], dupeTest[j]) && dupeTest[j].duplicate != true){
                //only true if it is both a duplicate and has not been previous marked as a duplicate
                dupeTest[j].duplicate = true;
                cout << i << ": " << j << endl;
                duplicates.push_back(dupeTest[i]);
                
            }else{
                //parents.push_back(dupeTest[i]);
            }
        }
    }
    //all that were not marked as duplicates are either unique or the original
    for(int i = 0; i < vectSize; i++){
        if(dupeTest[i].duplicate != true){
            parents.push_back(dupeTest[i]);
        }
    }
    //print the duplicates
    cout << "Duplicate size: " << duplicates.size() << endl;
    cout << "Duplicates index: " << endl;
    for(int i = 0; i < duplicates.size(); i++){
        cout << i << ": " << duplicates[i].posDiff << ", " << duplicates[i].speedDiff << endl;

    }
    //print the parents
    cout << "Parents size: " << parents.size() << endl;
    cout << "Parents index: " << endl;
    for(int i = 0; i < parents.size(); i++){
        cout << i << ": " << parents[i].posDiff << ", " << parents[i].speedDiff << endl;
    }
    //print the results
    for(int i = 0; i < dupeTest.size(); i++){
        if(dupeTest.size() > 50){//print every 50th adult for large vector
            if(i % 50 == 0){
                cout << "posDiff: " << dupeTest[i].posDiff << ", speedDiff: " << dupeTest[i].speedDiff << ", duplicate: " << dupeTest[i].duplicate << endl;
            }
        }
        else{
            cout << "posdiff: " << dupeTest[i].posDiff << ", speedDiff: " << dupeTest[i].speedDiff << ", duplicate: " << dupeTest[i].duplicate << endl;
        }
    }

    return true;

}