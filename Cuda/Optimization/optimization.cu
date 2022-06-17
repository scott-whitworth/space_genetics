// Didymos Optimization Project using CUDA and a genetic algorithm

//TODO: Clarify complexities of the include paths
//TODO: What / why we are including
#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
//#include "../Genetic_Algorithm/individuals.h" // For individual structs, paths to rkParameters for randomParameters()
#include "../Genetic_Algorithm/adult.h" // For adult structs, paths to rkParameters for randomParameters()
#include "../Genetic_Algorithm/child.h" // For child structs, paths to rkParameters for randomParameters()
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // for testing rk4simple
#include "../Genetic_Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration()
//
//#include "../Unit_Testing/testing_sorts.cpp"

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>   // allows us to use vectors instead of just arrays

//----------------------------------------------------------------------------------------------------------------------------
//Used to give rankings for sorting based on non-dominated sorting method. Used for rendezvous mission
//Assigns suitability rank to all adults.
//Input: pool - this generation of adults, defined/initilized in optimimize
//       cConstants
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Gives a distance value to each adult. A higher distance indicates that it is more diverse from other adults
// The distance is how different an adult is from the two adults closest to it for each objective function.
// Input: pool - the full pool of 2N adults excluding NaNs
//        poolSize - the poolSize of all the adults excluding NaNs
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Used to identify if the best adult has changed
// Input: previousBestPosDiff - the posDiff of the previous best adult
//        currentBest - the best Adult currently identified
//        distingushRate - utility variable used to help determine if the current best and previous best is different, will allow for comparison, regardless of how small the posDiffs are
// Output: Returns true if there was a change in the best individual
bool changeInBest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate);

//----------------------------------------------------------------------------------------------------------------------------
// ** Assumes pool is sorted array of Adults **
// Used in determining if main optimize loop continues
// Input: posTolerance - posDiff threshold, determines max target distance
//        speedTolerance - speedDiff threshold, determines max speed at target
//        pool - this generation of Adults, defined/initilized in optimimize
//        cConstants - struct holding config values, used for accessing best_count value
// Output: Returns true if top best_count adults within the pool are within the tolerance
bool allWithinTolerance(double posTolerance, double speedTolerance, const std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
//Function that will create the first generation of individuals within inputParamaeters
//Input: N set of adults as inputParameters, cConstants
//Output: new Adults of the first generation based on random parameters
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng);

//----------------------------------------------------------------------------------------------------------------------------
//Function that fill the allIndividuals array
//Input: allAdults, newAdults, oldAdults 
//Output: filled allAdults with the new and old
//TODO: Is the necessary?
//void fillAllAdults (std::vector<Adult> allAdults, std::vector<Adult> newAdults, std::vector<Adult> oldAdults, const cudaConstants* cConstants, const int& generation);

//----------------------------------------------------------------------------------------------------------------------------
//Function that will call the right sorting methods for the allIndividuals array
//Input: allAdults, numNans, cConstants 
//Output: impact sorted by posDiff, rendezvous sorted by rank and distance
void callSorts (std::vector<Adult>& allAdults, const int & numNans, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Function that will sort the adults from oldAdults and newAdults into allAdults
// Input: allAdults, oldAdults, and newAdults
// Output: oldAdults is filled with the adults that are potential parents for the next generation
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants, const int & generation);

//----------------------------------------------------------------------------------------------------------------------------
//Input: current anneal and dRate
//Output: an update of the anneal and dRate based on the tolerance and changeInBest
//Function that adjusts the anneal based on the current circumstances
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate);

//----------------------------------------------------------------------------------------------------------------------------
//This function will calculate the cost based on the mission type and the current best individual
//      This function should only be used within changeAnneal 
//Inputs: oldAdults - A vector of adults that will pull the individual whose final position will be used to calculate the cost
//              NOTE: This function assumes that oldAdults is already sorted by rankDistance and will pull the first adult from the vector
//Output: A cost double that will be used to change anneal              
double calculateCost(const std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// TEST / LIKELY TEMPORARY FUNCTION
// This function will find the minimum, maximum, and average distance of the individuals in allAdults, which will then be used for reporting
// 
// Inputs:  allAdults - array of adults that will be considered
//          minDist - minimum distance that will be calculated
//          avgDist - the average distance that will be calculated
//          maxDist - the maximum distance that will be calculated
// Outputs: The min, avg, and max distance variables will be filled in with the up-to-date values for this generation
void calculateDistValues (const std::vector<Adult> & allAdults, double & minDist, double & avgDist, double & maxDist);

//----------------------------------------------------------------------------------------------------------------------------
//Input: all the updated parameters of the current generation
//Output: calls the various terminal display functions when needed
//Function that handles the reporting of a generation's performance
void reportGeneration (std::vector<Adult>& oldAdults, const cudaConstants* cConstants, const double & currentAnneal, const double & anneal_min, const int & generation, int & numNans);

//----------------------------------------------------------------------------------------------------------------------------
// Main processing function for Genetic Algorithm
// - manages memory needs for genetic algorithm
// - deals with processing calls to CUDA callRK
// - exits when individuals converge on tolerance defined in Constants
double optimize(const cudaConstants* cConstants);


//Temp function that adds up the number of each status in an adults array and outputs the result
void countStatuses (const std::vector<Adult> & adultVec, const int & generation) {
    //Create ints for each status
    int numSuns, numNans, numValid, numOther;

    //Make them equal 0
    numSuns = numNans = numValid = numOther = 0;

    //Go thru the vector and add up the statuses
    for (int i = 0; i < adultVec.size(); i++)
    {
        if (adultVec[i].errorStatus == SUN_ERROR)
        {
            numSuns++;
        }
        else if (adultVec[i].errorStatus == NAN_ERROR)
        {
            numNans ++;
        }
        else if (adultVec[i].errorStatus == VALID)
        {
            numValid++;
        }
        else {
            numOther++;
        }
    }

    //Print the results
    std::cout << "\n\n_-_-_-_-_-_-_-_-_-_GENERATION #" << generation << " ERROR COUNTS_-_-_-_-_-_-_-_-_-_\n\n";

    std::cout << "\tSun Errors: " << numSuns;
    std::cout << "\n\tNan Errors: " << numNans;
    std::cout << "\n\tValids: " << numValid;
    std::cout << "\n\tOther: " << numOther;

    std::cout << "\n_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n";
    
}

//-----------------------------------------------------------------------------------------------------------------------------
int main () {
    
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "\n\nDevice Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    cudaSetDevice(0);

    // Declare the genetic constants used, with file path being used to receive initial values
    cudaConstants * cConstants = new cudaConstants("../Config_Constants/genetic.config"); 

    //test_main();

    // Sets run0 seed, used to change seed between runs
    // Seed is set in cudaConstants: current time or passed in via config
    double zero_seed = cConstants->time_seed;
    // Perform the optimization with optimize function
    for (int run = 0; run < cConstants->run_count; run++) {
        // Adjust the time_seed so it is unique based on each run
        cConstants->time_seed = zero_seed + run*100;

        // Display contents of cConstants being used for this run and how many runs
        std::cout << *cConstants;
        std::cout << "\tPerforming run #" << run+1 << "\n\n";

        // pre-calculate a table of Earth's position within possible mission time range
        // defined as global variable
        // accessed on the CPU when individuals are initilized
        launchCon = new EarthInfo(cConstants); 

        // File output of element values that were calculated in EarthInfo constructor for verification
        /*if (cConstants->record_mode == true) {
            recordEarthData(cConstants, run);
        }*/
        // Call optimize with the current parameters in cConstants
        optimize(cConstants);

        delete launchCon; // Deallocate launchCon info for this run as it may be using a different time range in the next run
    }
    // Now that the optimize function is done (assumed that optimize() also records it), deallocate memory of the cudaConstants
    delete cConstants;
    
    return 0;
}

//gives each adult in the allAdults vector a rank
void giveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants) {
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
    domination.resize(allAdults.size());

    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    std::vector<int> dominatedByCount;
    
    //fill the vector with 0s to make sure the count is accurate
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){

            //check the status of both i and j and see if i is automatically dominated
            if(allAdults[i].errorStatus != VALID && allAdults[j].errorStatus == VALID){
                dominatedByCount[i]++;

            }//check the status of both i and j and see if j is automatically dominated
            else if(allAdults[j].errorStatus != VALID && allAdults[i].errorStatus == VALID){
                domination[i].push_back(j);
                
            }//if either both are valid or both are not valid, it will rank them normally
            //Check to see if i dominates j
            else if (dominationCheck(allAdults[i], allAdults[j], cConstants)){
                //Put the jth index in the set of individuals dominated by i
                //std::cout << "\n" << i << "th (i) Adult dominates " << j << "th (j) Adult!\n";
                domination[i].push_back(j);
            }
            //Check to see if j dominates i
            else if ( dominationCheck( allAdults[j], allAdults[i], cConstants) ){
                //std::cout << "\n" << j << "th (j) Adult dominates " << i << "th (i) Adult!\n";
                //Add one to i's dominated by count
                dominatedByCount[i]++; 
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            front.push_back(i);

            //std::cout << "\n\nAdult #" << i << " ranked " << 1;
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

                    //std::cout << "\n\nAdult #" << domination[front[k]][l] << " ranked " << rankNum+1;

                    //Add the index of the this adult to newFront
                    newFront.push_back(domination[front[k]][l]);                        
                }
            }
        }
        //increment the rank number
        rankNum++;
        
        //empty the current front
        std::vector<int>().swap(front);

        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
    }
    //std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~FINISHED RANKING~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

//gives each adult in the allAdults a distance representing how different it is from other individuals
void giveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants){
    //pool that holds the indexes of the valid adults
    std::vector<int> validAdults;

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(allAdults.begin(), allAdults.begin() + allAdults.size(), rankSort);

    //checks if the adult is valid and then adds that index to the vector
    //the size of this vector will be used to find the distance for valid adults only
    for(int i = 0; i < allAdults.size(); i++){
        if(allAdults[i].errorStatus == VALID){
            validAdults.push_back(i);
        }
    }

    //set all to zero including the invalid adults
    for (int i = 0; i < allAdults.size(); i++ ){
        //reset each individual's distance
        allAdults[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(allAdults.begin(), allAdults.begin() + validAdults.size(), LowerPosDiff);
    //Set the boundaries
    allAdults[0].distance += MAX_DISTANCE; //+=1
    allAdults[validAdults.size() - 1].distance += MAX_DISTANCE; //+=1


    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    for(int i = 1; i < validAdults.size() - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffLeft = allAdults[i+1].posDiff/allAdults[validAdults.size() - 1].posDiff;
        normalPosDiffRight = allAdults[i-1].posDiff/allAdults[validAdults.size() - 1].posDiff;
        //distance = distance + abs((i+1) - (i-1))
        allAdults[i].distance = allAdults[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(allAdults[validAdults.size() - 1].posDiff - allAdults[0].posDiff));
    }

    //Repeat above process for speedDiff
    if(cConstants->missionType == Rendezvous){//only do this for the rendezvous mission since it has 2 objectives
        std::sort(allAdults.begin(), allAdults.begin() + validAdults.size(), LowerSpeedDiff);
        //Set the boundaries
        allAdults[0].distance += MAX_DISTANCE; //+=1
        allAdults[validAdults.size() - 1].distance += MAX_DISTANCE; //+=1

    
        //For each individual besides the upper and lower bounds, make their distance equal to
        //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
        double normalSpeedDiffLeft;
        double normalSpeedDiffRight;
        for(int i = 1; i < validAdults.size() - 1; i++){
            //Divide left and right individuals by the worst individual to normalize
            normalSpeedDiffLeft = allAdults[i+1].speedDiff/allAdults[validAdults.size() - 1].speedDiff;
            normalSpeedDiffRight = allAdults[i-1].speedDiff/allAdults[validAdults.size() - 1].speedDiff;
            //distance = distance + abs((i+1) - (i-1))
            allAdults[i].distance = allAdults[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(allAdults[validAdults.size() - 1].speedDiff - allAdults[0].speedDiff));
        }
    }
}

//tests if the best value has changed much in the since the last time this was called
bool changeInBest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate) {
    //truncate is used here to compare doubles via the distinguguishRate, to ensure that there has been relatively no change.
        if (trunc(previousBestPosDiff/distinguishRate) != trunc(currentBest.posDiff/distinguishRate)) {
            return true;
        }
        else { 
            return false;
        }
}

//Returns true if top best_count adults within the oldAdults vector are within the tolerance
bool allWithinTolerance(double posTolerance, double speedTolerance, const std::vector<Adult>& oldAdults, const cudaConstants* cConstants) {

    //Check what type of mission is running to use the correct posDiff function
    if (cConstants->missionType == Rendezvous){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            //both objectives need to be within tolerance
            //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CON POS CHECK-_-_-_-_-_-_-_-_-_\n\n";
            if (oldAdults[i].posDiff >= posTolerance){
                return false;
            }
            //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CON SPEED CHECK-_-_-_-_-_-_-_-_-_\n\n";
            if (oldAdults[i].speedDiff >= speedTolerance){
                return false;
            }
        }
    }
    else if(cConstants->missionType == Impact){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            
            if(oldAdults[i].posDiff >= posTolerance) {
                //One was not within 
                return false;
            }
        }  
    }

    // If iterated through and all were within tolerance, success
    return true;
}

//Creating a generation either randomly or based on values from a file
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng){
    Child* initialChildren = new Child[cConstants->num_individuals]; 
    // Initilize individuals randomly or from a file
    if (cConstants->random_start) {
        // individuals set to randomly generated, but reasonable, parameters
        for (int i = 0; i < cConstants->num_individuals; i++) { 
            initialChildren[i] = Child(randomParameters(rng, cConstants), cConstants);
        }
    }
    // Read from file using cConstants initial_start_file_address to get path
    else {
        // **Might be depreciated, not tested summer 2020**
        // Sets inputParameters to hold initial individuals based from file optimizedVector.bin
        const int numStarts = 14; // the number of different sets of starting parameters in the input file
        std::ifstream starts;
        starts.open(cConstants->initial_start_file_address, std::ifstream::in|std::ios::binary); // a file containing the final parameters of converged results from CPU calculations        

        // sort the data into 2 dimensions
        // one row is one set of starting parameters
        // each column is a specific variable:
        double startDoubles;
        // arrayCPU needs to be updated to handle the fact that OPTIM_VARS may be flexible
        double arrayCPU[numStarts][OPTIM_VARS];
        
        for (int i = 0; i < OPTIM_VARS; i++) { // rows
            for (int j = 0; j < numStarts; j++) { // columns
                starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
                arrayCPU[j][i] = startDoubles;
            }
        }
        starts.close();

         // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
        for (int i = 0; i < cConstants->num_individuals; i++) {
            int row = rng() % numStarts; // Choose a random row to get the parameters from

            double tripTime = arrayCPU[row][TRIPTIME_OFFSET];
            double alpha = arrayCPU[row][ALPHA_OFFSET];
            double beta = arrayCPU[row][BETA_OFFSET];
            double zeta = arrayCPU[row][ZETA_OFFSET];

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = arrayCPU[row][j + GAMMA_OFFSET];
            }

            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] =  arrayCPU[row][j + TAU_OFFSET];
            }

            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = arrayCPU[row][j + COAST_OFFSET];
            }

            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

            initialChildren[i] = Child(example, cConstants);
        }
    }
    firstGeneration(initialChildren, oldAdults, cConstants);
    delete[] initialChildren;
}

//Function that calls the sorting algorithims, depending on the selected mission
void callSorts (std::vector<Adult>&allAdults, const int & numNans, const cudaConstants* cConstants){
    /*if (cConstants->missionType == Impact) {
        //Decide the next generation of potential parents based on posDiff? was cost
        std::sort(allAdults.begin(), allAdults.end());
    }
    else if (cConstants->missionType == Rendezvous) {*/
        //give a rank to each adult based on domination sort
        //* Ignore any nans at the end of allAdults
        //must be called after checking for nans and before giveDistance
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE RANK-_-_-_-_-_-_-_-_-_\n\n";
        giveRank(allAdults, cConstants); //gives a rank to each adult
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE DISTANCE-_-_-_-_-_-_-_-_-_\n\n";
        giveDistance(allAdults, cConstants); //gives a distance to each adult
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE SORT-_-_-_-_-_-_-_-_-_\n\n";
        std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort); //sorts allAdults using rankDistance sort
    //} 
}

//fills oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants, const int & generation){
    std::vector<Adult>().swap(allAdults); //ensures this vector is empty and ready for new inputs
    numNans = 0;

    //TODO: Making sure there are no errors within allAdults is temporary, eventually we need to make sure they are at the end of the rankDistanceSort but still include them within allAdults

    //countStatuses(newAdults, generation);
    //countStatuses(oldAdults, generation);

    for (int i = 0; i < newAdults.size(); i++){ //copies all the elements of newAdults into allAdults
        if(newAdults[i].errorStatus != VALID){ //tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
        }
        else {
            allAdults.push_back(newAdults[i]);
        }
    }
    
    for (int i = 0; i < oldAdults.size(); i++){ //copies over all the elements of oldAdults into allAdults
        if(oldAdults[i].errorStatus != VALID){//tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
        }
        else {
            allAdults.push_back(oldAdults[i]);
        }
    }

    //countStatuses(allAdults, generation);

    //Sort the set of adults
    //callSorts(allAdults, numNans, cConstants);
    //give a rank to each adult based on domination sort
    //* Ignore any nans at the end of allAdults
    //must be called after checking for nans and before giveDistance
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE RANK-_-_-_-_-_-_-_-_-_\n\n";
    giveRank(allAdults, cConstants); //gives a rank to each adult
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE DISTANCE-_-_-_-_-_-_-_-_-_\n\n";
    giveDistance(allAdults, cConstants); //gives a distance to each adult
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE SORT-_-_-_-_-_-_-_-_-_\n\n";
    std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort); //sorts allAdults using rankDistance sort

    oldAdults.clear(); //empties oldAdults so new values can be put in it
    int counter = 0; //a variable that ensures we put the correct number of adults in oldAdults
    //copies the best adults from allAdults into oldAdults (should be half of allAdults that are copied over)
    //TODO:: Make sure dividing num_individuals / 2 is correct or if we should divide by something else

    while (counter < cConstants->num_individuals && counter < allAdults.size()){
        oldAdults.push_back(allAdults[counter]);
        counter++;
    }

    //error message prints if somehow there are not enough adults to fill oldIndividuals to the size it should be  
    if(counter == allAdults.size() && counter < cConstants->num_individuals){ //TODO: May or may not want this. If used before oldAdults and newAdults were both filled once, delete error message
        std::cout << "There are not enough adults to fill oldIndividuals properly" << std::endl;
    }

    //Clear the newAdult vector so that it is ready to recive new children
    newAdults.clear(); 

    //countStatuses(oldAdults, generation);
}

//TODO: Figure out what to replace posDiff (what used to be cost)
//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate){
    
    //Calculate the current cost for this generation
    double curCost = calculateCost(oldAdults, cConstants);
    
    //Caluclate the new current anneal
    //It will be a linear decrease from the initial/max anneal based on how well the current cost is compared to the cost threshold
    //If we happen to be close to covergence, and the cost we get is -0.0001:
  //currentAnneal =               0.1          * (              1.000001                  )
    //if the signs were correct, this would still be only a small decrease when very close to convergence
    currentAnneal = cConstants->anneal_initial * (1 - (cConstants->costThreshold / curCost));

    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
    

    /*
    // Scaling anneal based on proximity to tolerance
    // Far away: larger anneal scale, close: smaller anneal
    if (cConstants->missionType == Impact) {
        //Impact is only based on posDiff, so proximity-based annealing only relies on how close posDiff is to tolerance.
        new_anneal = currentAnneal * (1 - (posTolerance / oldAdults[0].posDiff));
    }

    else if (cConstants->missionType == Rendezvous) {
        if (posTolerance < oldAdults[0].posDiff){ 
            //TODO: decide what we want to do with this annealing   
            //Exponentially changing annealing TODO: this annealing and the impact annealing are hardly getting adjusted
            new_anneal = currentAnneal * (1 - pow(posTolerance / oldAdults[0].posDiff,2.0));
            if (new_anneal < cConstants->anneal_final){
                new_anneal = cConstants->anneal_final; //Set a true minimum for annealing
            }
        }
    }
    
    /*
    //Process to see if anneal needs to be adjusted
    // If generations are stale, anneal drops
    Adult currentBest;
    // Compare current best individual to that from CHANGE_CHECK (50) many generations ago.
    // If they are the same, change size of mutations
    if (static_cast<int>(generation) % (cConstants->change_check) == 0) { 
        currentBest = oldAdults[0];
        // checks for anneal to change
        // previousBest starts at 0 to ensure changeInBest = true on generation 0
        if ( !(changeInBest(previousBestPosDiff, currentBest, dRate)) ) { 
            //this ensures that changeInBest never compares two zeros, thus keeping dRate in relevance as the posDiff lowers
            if (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                while (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                    dRate = dRate/10; 
                }
                std::cout << "\nnew dRate: " << dRate << std::endl;
            }
            // If no change in BestIndividual across generations, reduce currentAnneal by anneal_factor while staying above anneal_min
            //reduce anneal_min
            anneal_min = cConstants->anneal_initial*exp(-sqrt(posTolerance/oldAdults[0].posDiff)*generation);
            if (anneal_min < cConstants->anneal_final){
                anneal_min = cConstants->anneal_final;//Set a true minimum for annealing
            }

            //Rendezvous mission uses anneal_min, Impact does not
            if(cConstants->missionType == Impact) {
                currentAnneal = currentAnneal * cConstants->anneal_factor;
            }
            else if (cConstants->missionType == Rendezvous){
                currentAnneal = (currentAnneal * cConstants->anneal_factor > anneal_min)? (currentAnneal * cConstants->anneal_factor):(anneal_min);
            }
            std::cout << "\nnew anneal: " << currentAnneal << std::endl;              
        }

        previousBestPosDiff = currentBest.posDiff; //posDiff here used to be cost
    }
    */
}

//Function that will calculate this generation's best adult's cost
double calculateCost(const std::vector<Adult> & oldAdults, const cudaConstants* cConstants){
    //Create the cost double that will be returned
    double cost; 

    //Variables will be used to calculate cost
    //      They both will specifically be used to calculate the adult's diff / goal diff
    //      This prevents a excellent speed diff causing the adult's diff / goal diff to be less than 1
    //          If this was the case, there is the potential that the cost goal could be met even if the position diff hadn't hit it's goal
    double rel_posDiff, rel_speedDiff;

    //How to calculate the cost will depend on the mission type
    //In general, cost will be (for each mission parameter, the difference/the goal) - the number of mission parameters
    //This will mean a cost of 0 or below signifies that the individual has hit the mission goals

    //For impacts, the cost only depends on posDiff
    //  NOTE: While the RD sorting favors high speed for impacts, convergence ultimately comes down to posDiff
    if (cConstants -> missionType == Impact) {
        //The goal is position threshold and the current status is the best individual's posDiff

        //Check to see if the rel_posDiff would be less than the threshold
        if (oldAdults[0].posDiff > cConstants->pos_threshold) {
            //If not, calculate how close the diff is to the goal
            rel_posDiff = oldAdults[0].posDiff / cConstants->pos_threshold; 
        }
        else {
            //If so, set the rel_posDiff to 1, signifying that the goal has been met
            rel_posDiff = 1.0;
        }

        //The coast is the rel_posDiff minus the number of mission goals (1 in this case)
        cost = rel_posDiff - IMPACT_MISSION_PARAMETER_COUNT;
    }
    //For rendezvous, the cost depends on both posDiff and speedDiff
    else {
        //Similarly to impact, calculate how far the best adult's speed and position diffs are away from the goal and subtract by the number of mission goals
        
        //First, same as above either set rel_posDiff to how close the best adult's posDiff is to the goal or to 1, depending on if the posDiff is beating the goal
        if (oldAdults[0].posDiff > cConstants->pos_threshold) {
            rel_posDiff = oldAdults[0].posDiff / cConstants->pos_threshold;
        }
        else {
            rel_posDiff = 1.0;
        }
        //Repeat the same process in calculating rel_posDiff for rel_speedDiff
        if (oldAdults[0].speedDiff > cConstants->speed_threshold) {
            rel_speedDiff = oldAdults[0].speedDiff / cConstants->speed_threshold;
        }
        else { 
            rel_speedDiff = 1.0;
        }

        //The cost the addition of each difference between the goals and pod/speed diffs minus the rendezvous mission goal count (2)
        //This means a flight that meets all goals will have a cost of 0, since the rel pos/speed diffs would be set to 1
        cost = (rel_posDiff + rel_speedDiff) - RENDEZVOUS_MISSION_PARAMETER_COUNT;
    }
    
    //Return the calculated cost
    return cost; 
}

//Function that will calculate distance values for a generation
void calculateDistValues (const std::vector<Adult> & allAdults, double & minDist, double & avgDist, double & maxDist){
    //Reset the dist values
    minDist = 2; //Set the min dist to the maximum possible value, so that it will be changed
    avgDist = 0; 
    maxDist = 0; //Set the max dist to the min possible value, so that it is garunteed to be changed
    
    //For loop will iterate through the adult array to find the values needed
    for (int i = 0; i < allAdults.size(); i++) {
        //Check to see if this adult's distance is a new minimum
        if (allAdults[i].distance < minDist) {
            minDist = allAdults[i].distance; //Set the new min distance
        }

        //Check to see if this adult's distance is the new maximum
        if (allAdults[i].distance > maxDist) {
            maxDist = allAdults[i].distance; //Ser the new max distance
        }

        
        
    }
    
}

//Function that writes the results of the inserted generation
void reportGeneration (std::vector<Adult> & oldAdults, const cudaConstants* cConstants, const double & currentAnneal, const double & anneal_min, const int & generation, int & numNans){
    // If in recording mode and write_freq reached, call the record method
    if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
        recordGenerationPerformance(cConstants, oldAdults, generation, currentAnneal, cConstants->num_individuals, anneal_min);
    }

    // Only call terminalDisplay every DISP_FREQ, not every single generation
    if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
        // Prints the best individual's posDiff / speedDiff

        //best position individual
        std::cout << "\n\nBest Position Individual:";
        std::sort(oldAdults.begin(), oldAdults.end(), LowerPosDiff);
        terminalDisplay(oldAdults[0], generation);

        if(cConstants->missionType == Rendezvous){
            //Best lower speed individual
            std::cout << "\nBest Speed Individual:";
            std::sort(oldAdults.begin(), oldAdults.end(), LowerSpeedDiff);
            terminalDisplay(oldAdults[0], generation);
        }
        else if(cConstants->missionType == Impact){
            //Best higher speed individual
            std::cout << "\nBest Speed Individual:";
            std::sort(oldAdults.begin(), oldAdults.end(), HigherSpeedDiff);
            terminalDisplay(oldAdults[0], generation);
        }
        //display to the terminal the best individual based on cost
        std::cout << "\nBest Rank Distance Individual:";
        std::sort(oldAdults.begin(), oldAdults.end(), rankDistanceSort);
        terminalDisplay(oldAdults[0], generation);
        std::cout << "\n# of errors this generation: " << numNans << "\n" << std::endl;
        
        //Reset the tally of nans.
        numNans = 0;
    }

    //Record the parent pool for the next generation
    if (static_cast<int>(generation) % cConstants->all_write_freq == 0 && cConstants->record_mode == true) {
        recordAllIndividuals("NextParents", cConstants, oldAdults, cConstants->num_individuals, generation);
    }
}

//Function that will facilitate the process of finding an optimal flight path
double optimize(const cudaConstants* cConstants) {
    // Not used, previously used for reporting computational performance
    double calcPerS = 0;

    time_t timeSeed = cConstants->time_seed;
    std::mt19937_64 rng(timeSeed); // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;

    // Initialize the recording files if in record mode
    if (cConstants->record_mode == true) {
        initializeRecord(cConstants);
    }

    // input parameters for Runge Kutta process
    // Each parameter is the same for each thread on the GPU
    double timeInitial = 0; // the starting time of the trip is always defined as zero   

    // Runge Kutta adaptive time step error tolerance
    // Not used now that callRK has been moved to ga_crossover
    //double absTol = cConstants->rk_tol; 

    // the starting step size for RK run
    // - note that the current step size varies throughout each run
    // - the inital stepSize greatly influences the number of steps the RK calculation will take
    // - previously, in place of max_numsteps we used a very large divisor GuessMaxPossibleSteps (resulting in very small step size) 
    //   this resulted in most, if not all of the rk calculations taking max_numsteps loops always
    double stepSize = ((cConstants->orbitalPeriod) - timeInitial) / cConstants->max_numsteps; 

    // Initial genetic anneal scalar
    double currentAnneal = cConstants->anneal_initial;

    //lower bound for anneal so it does not get too small. Only used with rendezvous mission.
    double anneal_min = cConstants->anneal_initial;
    
    // Main set of parameters for Genetic Algorithm
    // contains all thread unique input parameters
    // The "children pool" of the current genertation
    std::vector<Adult> newAdults; 

    //Input parameters from the previous generation. Mixes with new children to determine new inputParameters
    // The "potential parent pool" of the current generation
    // DISCLAIMER - this is mentioned as the "parent pool", but it's actually the *POTENTIAL* parent pool. The survivor pool is what's used to generate the next generation. Survivors are taken from this, so it's more accurate to call this the "Potential Parent Pool"
    std::vector<Adult> oldAdults; 

    //the set of all old and new individuals
    std::vector<Adult> allAdults;

    //used to compare the progress in cost between generations
    //Set as 0 at first to make sure that there is a change in the best indvidivual for the first check
    double previousBestPosDiff = 0;

    // number of current generation
    double generation = 0;    
    
    // Genetic solution tolerance for each objective
    double posTolerance = cConstants->pos_threshold;
    double speedTolerance = cConstants->speed_threshold;  
             
    // distinguishable rate used in changeInBest()
    //  - used to help check for a change in anneal
    //  - Gets smaller when no change is detected
    double dRate = 1.0e-8;

    // Flag for finishing the genetic process
    // set by allWithinTolerance()
    bool convergence = false;

    //sets the anneal for the generation
    //used for proximity based annealing
    double new_anneal;

    //A place to hold the number of Nans for recording at the end - numNans is for all adults in the generation - the oldAdults and the newAdults that have nan values are both tallied here
    int numNans = 0;

    //Creates the individuals needed for the 0th generation
    //Need to make children, then callRK, then make into adults (not currently doing that)
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE FIRST GEN-_-_-_-_-_-_-_-_-_\n\n";
    createFirstGeneration(oldAdults, cConstants, rng); 

    // main gentic algorithm loop
    // - continues until allWithinTolerance returns true (specific number of individuals are within threshold)
    do {        
        // Genetic Crossover and mutation occur here
        //takes in oldAdults (the potential parents) and fills newAdults with descendants of the old adults
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE NEW GEN-_-_-_-_-_-_-_-_-_\n\n";
        newGeneration(oldAdults, newAdults, currentAnneal, generation, rng, cConstants);

        /*
        TEST 
        for (int i = 0; i < 10; i++)
        {
            int randIndex = rng() % (cConstants -> num_individuals);

            std::cout << "Rand Adult (adult #" << randIndex << ") pos and speed diffs: \n";
            std::cout << "\tposDiff: " << newAdults[i].posDiff << "\n";
            std::cout << "\tspeedDiff: " << newAdults[i].speedDiff << "\n\n";
        }
        */
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE PREP PARENTS-_-_-_-_-_-_-_-_-_\n\n";
        //fill oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected (numNans is for all adults in the generation - the oldAdults and the newAdults)
        preparePotentialParents(allAdults, newAdults, oldAdults, numNans, cConstants, generation);

        //TODO:: Major space for efficeincy change
        /*
            in reference to 2002-Deb_NSGA-II.pdf, 186/page 5:
            old individuals (Qt) are added to the newly computed next generation (Pt)
            The issue is that Qt has already been computed, we should look to see if we can limit the amount we re-calculate individuals
            This would mean we need to store if it has been calculated or not

            Possible direction: just put callRK inside newGeneration() (ga_crossover.cpp)
        */
        //UPDATE (6/7/22): callRK moved into ga_crossover, the effect on the efficiency of the code is not known yet
        //Output survivors for the current generation if write_freq is reached
        if (static_cast<int>(generation) % cConstants->all_write_freq == 0 && cConstants->record_mode == true) {
            recordAllIndividuals("Survivors", cConstants, oldAdults, cConstants->survivor_count, generation);
        }

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the terminalDisplay() calls across generations 
        std::cout << '.';

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE ANNEAL STFF-_-_-_-_-_-_-_-_-_\n\n";
        //Perform utitlity tasks (adjusting anneal and reporting data)
        //Assumes oldAdults is in some order
        changeAnneal (oldAdults, cConstants, new_anneal, currentAnneal, anneal_min, previousBestPosDiff, generation, posTolerance, dRate);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE RECORD-_-_-_-_-_-_-_-_-_\n\n";
        reportGeneration (oldAdults, cConstants, currentAnneal, anneal_min, generation, numNans);

        // Before replacing new adults, determine whether all are within tolerance
        // Determines when loop is finished
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CONVERGENCE CHECK-_-_-_-_-_-_-_-_-_\n\n";
        convergence = allWithinTolerance(posTolerance, speedTolerance, oldAdults, cConstants);
        
        //Increment the generation counter
        ++generation;

        //sorts oldAdults using rankDistance sort
        //This will prep the oldAdults array to be used to generate new children; this makes sure any sorts from the record/report functions are undone
        std::sort(oldAdults.begin(), oldAdults.end(), rankDistanceSort); 
    
        //Loop exits based on result of allWithinTolerance and if max_generations has been hit
    } while ( !convergence && generation < cConstants->max_generations);

    // Call record for final generation regardless of frequency
    // for the annealing argument, set to -1 (since the anneal is only relevant to the next generation and so means nothing for the last one)
    // for the numFront argument, set to -1 (just because)
    if (cConstants->record_mode == true) {
        recordGenerationPerformance(cConstants, oldAdults, generation, currentAnneal, cConstants->num_individuals, anneal_min);
    }
    // Only call finalRecord if the results actually converged on a solution
    // also display last generation onto terminal
    if (convergence) {
        terminalDisplay(oldAdults[0], generation);
        finalRecord(cConstants, oldAdults, static_cast<int>(generation));
    }

    return calcPerS;
}