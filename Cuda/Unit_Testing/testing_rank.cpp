#include "../Genetic_Algorithm/adult.h"
#include "../Genetic_Algorithm/individuals.h"
//#include UNITTEST
#include <iostream>
#include <iomanip>
#include <vector>
#include <string> //allows us to use to_string and string
using std::cout;
using std::endl;


// this function will check if giveRank is functioning as expected
bool dominationTest();

//test version of the giveRank function from optimization
void giveRankTest(std::vector<Adult> & allAdults);

//test version of the giveDistance function from optimization
void giveDistanceTest(std::vector<Adult> pool, int poolSize);

//giveRank fucntion that worked using arrays from individual
void oldGiveRankTest(Individual * pool, int poolSize);

//for testing the changeAnneal and changeInBest functions
bool CAtest();

//test version of changeAnneal from optimization.cu
void changeAnnealTest(std::vector<Adult> newAdults, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate);

bool changeInBestTest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate);

bool dominationTest(){

    //vector for the test
    std::vector<Adult> GRtest;

    int r = 0;//rank
    double d = 1.0;//distance
    int vectSize = 6;//size of the vector for the test

    //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        GRtest.push_back(Adult(r, d)); 
        //GRtest[i].posDiff = i;
        //GRtest[i].speedDiff = i;
    }
    GRtest[0].posDiff = 2;
    GRtest[0].speedDiff = 4;
    //should be rank 4; expected dominatedByCount = 2
    GRtest[1].posDiff = 4;
    GRtest[1].speedDiff = 1;
    //should be rank 3; expected dominatedByCount = 1?
    GRtest[2].posDiff = 5;
    GRtest[2].speedDiff = 7;
    //should be rank 2; expected dominatedByCount = 1
    GRtest[3].posDiff = 9.72;
    GRtest[3].speedDiff = 9.72;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[4].posDiff = 11;
    GRtest[4].speedDiff = 2;
    //should be rank 1; expected dominatedByCount = 0
    GRtest[5].posDiff = 3;
    GRtest[5].speedDiff = 3;
    //
    //GRtest[6].posDiff = 10;
    //GRtest[6].speedDiff = 10;
    //
    //GRtest[7].posDiff = 5;
    //GRtest[7].speedDiff = 5;
    
    
    giveRankTest(GRtest);
    cout << "Ranks:" << endl;
    for(int i = 0; i < GRtest.size(); i++){
        if(i % 1 == 0){
           cout << "r: " << GRtest[i].rank << ", p: " << GRtest[i].posDiff << ", s: " << GRtest[i].speedDiff << endl; 
        }
          
    }
    
    giveDistanceTest(GRtest, GRtest.size());


    Individual * testPool = new Individual[vectSize];

        //creating the vector and filling it with adults
    for(int i = 0; i < vectSize; i++){
        testPool[i] = Individual(r, d); 
        testPool[i].posDiff = i;
        testPool[i].speedDiff = i;
    }
    /*
    //giving each adult specific posDiffs and speedDiffs for expected results
    //should be rank 5; expected dominatedByCount = 4
    testPool[0].posDiff = 25;
    testPool[0].speedDiff = 25;
    //should be rank 4; expected dominatedByCount = 2
    testPool[1].posDiff = 15;
    testPool[1].speedDiff = 15;
    //should be rank 3; expected dominatedByCount = 1?
    testPool[2].posDiff = 1;
    testPool[2].speedDiff = 1;
    //should be rank 2; expected dominatedByCount = 1
    testPool[3].posDiff = 1;
    testPool[3].speedDiff = 1;
    //should be rank 1; expected dominatedByCount = 0
    testPool[4].posDiff = 1;
    testPool[4].speedDiff = 1;
    //should be rank 1; expected dominatedByCount = 0
    testPool[5].posDiff = 1;
    testPool[5].speedDiff = 1;
    //
    testPool[6].posDiff = 25;
    testPool[6].speedDiff = 25;
    //
    testPool[7].posDiff = 25;
    testPool[7].speedDiff = 25;
    

    oldGiveRankTest(testPool, vectSize);

    cout << "Ranks:" << endl;
    for(int i = 0; i < vectSize; i++){
        cout << "r: " << testPool[i].rank << ", p: " << testPool[i].posDiff << ", s: " << testPool[i].speedDiff << endl;  
    }
    
    */
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
    domination.resize(allAdults.size(), std::vector<int>(allAdults.size()));
    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    //Note: fill the vector with 0s to make sure the count is accurate
    //TODO: unit test to make sure the whole vector is actually initially filled with 0's and not just the first index or the original vector size
    std::vector<int> dominatedByCount; 
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){
        
        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){

            //Check to see if i dominates j
            if (dominationCheckTest(allAdults[i], allAdults[j])){
                //Put the jth index in the set of individuals dominated by i
                domination[i][j] = j;
            }
            //Check to see if j dominates i
            else if (dominationCheckTest(allAdults[j], allAdults[i])) {
                dominatedByCount[i]++; 
                
            }
            
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            
            //cout << " front size: " << front.size() << endl;
            front.push_back(i);
        }
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

    for(int i = 0; i < dominatedByCount.size(); i++){
            cout << dominatedByCount[i] << ", ";
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
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[front[k]][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[front[k]][l]].rank = rankNum + 1;
                    //Add the index of the this adult to newFront

                    newFront.push_back(l);
                }
                
            }
        }

        //increment the rank number
        rankNum++;
        //cout << " current overall rank: " << rankNum << endl;
        //empty the current front
        std::vector<int>().swap(front);

        cout << "new counts: " << endl;
        for(int i = 0; i < dominatedByCount.size(); i++){
            cout << dominatedByCount[i] << ", ";
        }
        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        cout << endl << "new front: " << endl;
        for(int i = 0; i < newFront.size(); i++){
            cout << newFront[i] << ", ";
        }
        front = newFront;
        cout << endl;
    }
}


//gives each adult in the pool a distance representing how different it is from other individuals
void giveDistanceTest(std::vector<Adult> pool, int poolSize){

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(pool.begin(), pool.begin() + poolSize, rankSort);

    cout << endl << "Post rank sort:" << endl;
    for(int i = 0; i < pool.size(); i++){
        if(i % 1 == 0){
           cout << "r: " << pool[i].rank << ", p: " << pool[i].posDiff << ", s: " << pool[i].speedDiff << endl; 
        }
          
    }
    for (int i = 0; i < poolSize; i++ ){
        //reset each individual's distance
        pool[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(pool.begin(), pool.begin() + poolSize, LowerPosDiff);
    //Set the boundaries
    //these dont make much sense, this makes the best posDiff in rank 1 be the worst overall in rank 1
    //pool[0].distance = 1e+12;
    //pool[poolSize - 1].distance = 1e+12;
    cout << endl << "Post lowerPosDiff sort:" << endl;
    for(int i = 0; i < pool.size(); i++){
        if(i % 1 == 0){
           cout << "r: " << pool[i].rank << ", p: " << pool[i].posDiff << ", s: " << pool[i].speedDiff << " d:" << endl;
        }
          
    }

    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    //could we instead have it find the distance from i = 0 and i = 1
    for(int i = 0; i < poolSize; i++){
        //Divide left and right individuals by the worst individual to normalize
        if(i == 0){
            normalPosDiffLeft = pool[i+1].posDiff/pool[poolSize - 1].posDiff;
            normalPosDiffRight = 0;
        }else if(i == poolSize - 1){
            normalPosDiffLeft = 1e+12;
            normalPosDiffRight = pool[i-1].posDiff/pool[poolSize - 1].posDiff;
        }else{
            normalPosDiffRight = pool[i-1].posDiff/pool[poolSize - 1].posDiff;
            normalPosDiffLeft = pool[i+1].posDiff/pool[poolSize - 1].posDiff;
        }
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(pool[poolSize - 1].posDiff - pool[0].posDiff));
        //cout << std::setprecision(4) << p << endl;
        cout << i << " posDiff Distance " <<  pool[i].distance << endl;
    }

    //Repeat above process for speedDiff    
    std::sort(pool.begin(), pool.begin() + poolSize, LowerSpeedDiff);
    //Set the boundaries
    //pool[0].distance = 1e+12;
    //pool[poolSize - 1].distance = 1e+12;
    cout << endl << "Post lowerSpeedDiff sort:" << endl;
    for(int i = 0; i < pool.size(); i++){
        if(i % 1 == 0){
           cout << "r: " << pool[i].rank << ", p: " << pool[i].posDiff << ", s: " << pool[i].speedDiff << " d:" << pool[i].distance << endl; 
        }
          
    }

    
    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalSpeedDiffLeft;
    double normalSpeedDiffRight;
    for(int i = 0; i < poolSize; i++){
        //Divide left and right individuals by the worst individual to normalize
        if(i == 0){
            normalSpeedDiffLeft = pool[i+1].speedDiff/pool[poolSize - 1].speedDiff;
            normalSpeedDiffLeft = 0;
        }else if(i == poolSize - 1){
            normalSpeedDiffLeft = 1e+12;
            normalSpeedDiffRight = pool[i-1].speedDiff/pool[poolSize - 1].speedDiff;
        }else{
            normalSpeedDiffRight = pool[i-1].speedDiff/pool[poolSize - 1].speedDiff;
            normalSpeedDiffLeft = pool[i+1].speedDiff/pool[poolSize - 1].speedDiff;
        }
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(pool[poolSize - 1].speedDiff - pool[0].speedDiff));

    }

    cout << endl << " total Distance " << endl;
    for(int i = 0; i < pool.size(); i++){
        if(i % 1 == 0){
           cout << "r: " << pool[i].rank << ", p: " << pool[i].posDiff << ", s: " << pool[i].speedDiff << " d:" << pool[i].distance << endl; 
        }
          
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
            if (dominationCheckTestv2(pool[i], pool[j])){
                pool[i].dominated.push_back(j);
            }
            //if j dominates i, increase the number of times that i has been dominated
            else if (dominationCheckTestv2(pool[j], pool[i])) {
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
    CAtest[0].posDiff = 1e-6;
    CAtest[0].speedDiff = 1;
    //should be rank 4; expected dominatedByCount = 2
    CAtest[1].posDiff = 4;
    CAtest[1].speedDiff = 1;
    //should be rank 3; expected dominatedByCount = 1?
    CAtest[2].posDiff = 5;
    CAtest[2].speedDiff = 7;
    //should be rank 2; expected dominatedByCount = 1
    CAtest[3].posDiff = 9.72;
    CAtest[3].speedDiff = 9.72;
    //should be rank 1; expected dominatedByCount = 0
    CAtest[4].posDiff = 11;
    CAtest[4].speedDiff = 2;
    //should be rank 1; expected dominatedByCount = 0
    CAtest[5].posDiff = 3;
    CAtest[5].speedDiff = 3;

    giveRankTest(CAtest);

    giveDistanceTest(CAtest, CAtest.size());

    double anneal_initial = 0.10;
    double new_anneal;
    double currentAnneal = anneal_initial;
    double anneal_min = anneal_initial;
    double previousBestPosDiff = 0;
    double generation = 1;
    double posTolerance = 1e-14;
    double dRate = 1e-8;

    while(generation < 10){
        changeAnnealTest(CAtest, new_anneal, currentAnneal, anneal_min, previousBestPosDiff, generation, posTolerance, dRate);
        cout << "done with run " << generation << endl;
        cout << " new_ anneal: " << new_anneal << " currentAnneal: " << currentAnneal << " anneal_min: " << anneal_min << endl;
        cout <<  " previousBestPosDiff: " << previousBestPosDiff << " dRate: " << dRate <<  endl;
        CAtest[0].posDiff = CAtest[0].posDiff/10;
        generation++;
    }

    return true;

}


void changeAnnealTest(std::vector<Adult> newAdults, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate){
    // Scaling anneal based on proximity to tolerance
    //cudaConstants:
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