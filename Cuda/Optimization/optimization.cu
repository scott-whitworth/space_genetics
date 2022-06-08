// Didymos Optimization Project using CUDA and a genetic algorithm

//TODO: Clarify complexities of the include paths
//TODO: What / why we are including
#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/individuals.h" // For individual structs, paths to rkParameters for randomParameters()
#include "../Genetic_Algorithm/adults.h" // For adult structs, paths to rkParameters for randomParameters()
#include "../Genetic_Algorithm/children.h" // For child structs, paths to rkParameters for randomParameters()
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // for testing rk4simple
#include "../Genetic_Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration()

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>   // allows us to use vectors instead of just arrays

//----------------------------------------------------------------------------------------------------------------------------
//Used to give rankings for sorting based on non-dominated sorting method. Used for rendezvous mission
//Assigns suitability rank to all adults.
//Input: pool - this generation of adults, defined/initilized in optimimize
//       cConstants
void giveRank(std::vector<Adult> pool, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Gives a distance value to each adult. A higher distance indicates that it is more diverse from other adults
// The distance is how different an adult is from the two adults closest to it for each objective function.
// Input: pool - the full pool of 2N adults excluding NaNs
//        poolSize - the poolSize of all the adults excluding NaNs
void giveDistance(std::vector<Adult> pool, const cudaConstants* cConstants, int poolSize);

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
bool allWithinTolerance(double posTolerance, double speedTolerance, std::vector<Adult> pool, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
//Function that will create the first generation of individuals within inputParamaeters
//Input: N set of adults as inputParameters, cConstants
//Output: new Adults of the first generation based on random parameters
void createFirstGeneration(std::vector<Adult> oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng);

//----------------------------------------------------------------------------------------------------------------------------
//Function that fill the allIndividuals array
//Input: allAdults, newAdults, oldAdults 
//Output: filled allAdults with the new and old
void fillAllAdults (std::vector<Adult> allAdults, std::vector<Adult> newAdults, std::vector<Adult> oldAdults, const cudaConstants* cConstants, const int& generation);

//----------------------------------------------------------------------------------------------------------------------------
//Function that will call the right sorting methods for the allIndividuals array
//Input: allAdults, numNans, cConstants 
//Output: impact sorted by posDiff, rendezvous sorted by rank and distance
void callSorts (std::vector<Adult> allAdults, const int & numNans, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Function that will sort the adults from oldAdults and newAdults into allAdults
// Input: allAdults, oldAdults, and newAdults
// Output: oldAdults is filled with the adults that are potential parents for the next generation
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
//Input: current anneal and dRate
//Output: an update of the anneal and dRate based on the tolerance and changeInBest
//Function that adjusts the anneal based on the current circumstances
void changeAnneal (std::vector<Adult> newAdults, const cudaConstants* cConstants, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate);

//----------------------------------------------------------------------------------------------------------------------------
//Input: all the updated parameters of the current generation
//Output: calls the various terminal display functions when needed
//Function that handles the reporting of a generation's performance
void reportGeneration (std::vector<Adult> newAdults, const cudaConstants* cConstants, const double & new_anneal, const double & anneal_min, const int & generation, int & numNans);

//----------------------------------------------------------------------------------------------------------------------------
// Main processing function for Genetic Algorithm
// - manages memory needs for genetic algorithm
// - deals with processing calls to CUDA callRK
// - exits when individuals converge on tolerance defined in Constants
double optimize(const cudaConstants* cConstants);

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

//TODO: unit test to make sure every individual is given NEW rank

//gives each adult in the pool a rank
void giveRank(std::vector<Adult> pool, const cudaConstants* cConstants) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

    //Used to store the current front of individuals. first filled with the first front individuals(best out of all population)
    // filled with index of individuals in pool
    std::vector<int> front;

    //TODO: Pull Adult::dominates and Adult::dominatedByCount into this function
    // probably a 2D vector, or an array of vectors
    
    //loop through each individual
    for (int i = 0; i < pool.size(); i++){
        
        //number of times pool[i] has been dominated
        pool[i].dominatedByCount = 0;

        //set of solutions that pool[i] dominates. Need to empty for each generation
        std::vector<int>().swap(pool[i].dominates);

        for(int j = 0; j < pool.size(); j++){
            
            //if i dominates j, put the j index in the set of individuals dominated by i.
            if (dominates(pool[i], pool[j], cConstants)){
                pool[i].dominates.push_back(j);
            }
            //if j dominates i, increase the number of times that i has been dominated
            else if (dominates(pool[j], pool[i], cConstants)) {
                pool[i].dominatedByCount++;
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (pool[i].dominatedByCount == 0){
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
            for(int l = 0; l < pool[front[k]].dominates.size(); l++){

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                pool[pool[front[k]].dominates[l]].dominatedByCount--;

                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (pool[pool[front[k]].dominates[l]].dominatedByCount == 0){
                    pool[pool[front[k]].dominates[l]].rank = rankNum + 1;
                    newFront.push_back(pool[front[k]].dominates[l]);                        
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

//gives each adult in the pool a distance representing how different it is from other individuals
void giveDistance(std::vector<Adult> pool, const cudaConstants* cConstants, int poolSize){

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(pool.begin(), pool.begin() + cConstants->num_individuals*2, rankSort);

    for (int i = 0; i < poolSize; i++ ){
        //reset each individual's distance
        pool[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(pool.begin(), pool.begin() + poolSize, LowerPosDiff);
    //Set the boundaries
    pool[0].distance = cConstants->MAX_DISTANCE;
    pool[poolSize - 1].distance = cConstants->MAX_DISTANCE;

    //TODO:: Constant set for these numbers MAX_DISTANCE

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
    }

    //Repeat above process for speedDiff    
    std::sort(pool.begin(), pool.begin() + poolSize, LowerSpeedDiff);
    //Set the boundaries
    pool[0].distance = cConstants->MAX_DISTANCE;
    pool[poolSize - 1].distance = cConstants->MAX_DISTANCE;

    
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

//Returns true if top best_count adults within the pool are within the tolerance
bool allWithinTolerance(double posTolerance, double speedTolerance, std::vector<Adult>pool, const cudaConstants* cConstants) {

    //Check what type of mission is running to use the correct posDiff function
    if (cConstants->missionType == Rendezvous){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            //both objectives need to be within tolerance
            if (pool[i].posDiff >= posTolerance){
                return false;
            }
            if (pool[i].speedDiff >= speedTolerance){
                return false;
            }
        }
    }
    else if(cConstants->missionType == Impact){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            
            if(pool[i].posDiff >= posTolerance) {
                //One was not within 
                return false;
            }
        }  
    }

    // If iterated through and all were within tolerance, success
    return true;
}

//Creating a generation either randomly or based on values from a file
void createFirstGeneration(std::vector<Adult>oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng){
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
void callSorts (std::vector<Adult>allAdults, const int & numNans, const cudaConstants* cConstants){
    if (cConstants->missionType == Impact) {
        //Decide the next generation of potential parents based on posDiff? was cost
        std::sort(allAdults.begin(), allAdults.end());
    }
    else if (cConstants->missionType == Rendezvous) {
        //give a rank to each adult based on domination sort
        //* Ignore any nans at the end of allAdults
        //must be called after checking for nans and before giveDistance
        giveRank(allAdults, cConstants); //gives a rank to each adult
        giveDistance(allAdults, cConstants, allAdults.size() - numNans); //gives a distance to each adult
        std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort); //sorts allAdults using rankDistance sort
    } 
}

//fills oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, const cudaConstants* cConstants){
    allAdults.clear(); //ensures this vector is empty and ready for new inputs
    numNans = 0;

    for (int i = 0; i < newAdults.size(); i++){ //copies all the elements of newAdults into allAdults
        allAdults.push_back(newAdults[i]);
        if(newAdults[i].status != VALID){ //tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
        }
    }
    for (int i = 0; i < oldAdults.size(); i++){ //copies over all the elements of oldAdults into allAdults
        allAdults.push_back(oldAdults[i]);
        if(oldAdults[i].status != VALID){//tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
        }
    }

    //Sort the set of adults
    callSorts(allAdults, numNans, cConstants);

    oldAdults.clear(); //empties oldAdults so new values can be put in it
    int counter = 0; //a variable that ensures we put the correct number of adults in oldAdults
    //copies the best adults from allAdults into oldAdults (should be half of allAdults that are copied over)
    //TODO:: Make sure dividing num_individuals / 2 is correct or if we should divide by something else
    while (counter < (cConstants->num_individuals)/2 && counter < allAdults.size()){
        oldAdults.push_back(allAdults[counter]);
        counter++;
    }

    //error message prints if somehow there are not enough adults to fill oldIndividuals to the size it should be  
    if(counter == allAdults.size() && counter < (cConstants->num_individuals)/2){ //TODO: May or may not want this. If used before oldAdults and newAdults were both filled once, delete error message
        std::cout << "There are not enough adults to fill a oldIndividuals properly" << std::endl;
    }

}

//TODO: Figure out what to replace posDiff (what used to be cost)
//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (std::vector<Adult> newAdults, const cudaConstants* cConstants, double & new_anneal, double & currentAnneal, double & anneal_min,  double & previousBestPosDiff, double & generation, const double & posTolerance, double & dRate){
    // Scaling anneal based on proximity to tolerance
    // Far away: larger anneal scale, close: smaller anneal
    if (cConstants->missionType == Impact) {
        //Impact is only based on posDiff, so proximity-based annealing only relies on how close posDiff is to tolerance.
        new_anneal = currentAnneal * (1 - (posTolerance / newAdults[0].posDiff));
    }

    else if (cConstants->missionType == Rendezvous) {
        if (posTolerance < newAdults[0].posDiff){ 
            //TO DO: decide what we want to do with this annealing   
            //Exponentially changing annealing, as oppose to what?
            new_anneal = currentAnneal * (1 - pow(posTolerance / newAdults[0].posDiff,2.0));
            if (new_anneal < cConstants->anneal_final){
                new_anneal = cConstants->anneal_final; //Set a true minimum for annealing
            }
        }
    }
    
    //Process to see if anneal needs to be adjusted
    // If generations are stale, anneal drops
    Adult currentBest;
    // Compare current best individual to that from CHANGE_CHECK (50) many generations ago.
    // If they are the same, change size of mutations
    if (static_cast<int>(generation) % (cConstants->change_check) == 0) { 
        currentBest = newAdults[0];
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
            anneal_min = cConstants->anneal_initial*exp(-sqrt(posTolerance/newAdults[0].posDiff)*generation);
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
}

//Function that writes the results of the inserted generation
void reportGeneration (std::vector<Adult> newAdults, const cudaConstants* cConstants, const double & new_anneal, const double & anneal_min, const int & generation, int & numNans){
    // If in recording mode and write_freq reached, call the record method
    if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
        recordGenerationPerformance(cConstants, newAdults, generation, new_anneal, cConstants->num_individuals, anneal_min);
    }

    // Only call terminalDisplay every DISP_FREQ, not every single generation
    if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
        // Prints the best individual's posDiff / speedDiff

        //best position individual
        std::cout << "\n\nBest Position Individual:";
        std::sort(newAdults.begin(), newAdults.end(), LowerPosDiff);
        terminalDisplay(newAdults[0], generation);

        if(cConstants->missionType == Rendezvous){
            //Best lower speed individual
            std::cout << "\nBest Speed Individual:";
            std::sort(newAdults.begin(), newAdults.end(), LowerSpeedDiff);
            terminalDisplay(newAdults[0], generation);
        }
        else if(cConstants->missionType == Impact){
            //Best higher speed individual
            std::cout << "\nBest Speed Individual:";
            std::sort(newAdults.begin(), newAdults.end(), HigherSpeedDiff);
            terminalDisplay(newAdults[0], generation);
        }
        //display to the terminal the best individual based on cost
        std::cout << "\nBest PosDiff Individual:";
        std::sort(newAdults.begin(), newAdults.end(), rankDistanceSort);
        terminalDisplay(newAdults[0], generation);
        std::cout << "\n# of Nans this generation: " << numNans << "\n" << std::endl;
        
        //re-sort by rankDistance for rendezvous mission
        if(cConstants->missionType == Rendezvous) {
            std::sort(newAdults.begin(), newAdults.end(), rankDistanceSort);
        }
        
        //Reset the tally of nans.
        numNans = 0;
    }

    //Record the parent pool for the next generation
    if (static_cast<int>(generation) % cConstants->all_write_freq == 0 && cConstants->record_mode == true) {
        recordAllIndividuals("NextParents", cConstants, newAdults, cConstants->num_individuals, generation);
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
    //TODO: Should this be based on max_numsteps?
    double stepSize = ((cConstants->orbitalPeriod) - timeInitial) / cConstants->GuessMaxPossibleSteps; 

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
    createFirstGeneration(oldAdults, cConstants, rng); 

    // main gentic algorithm loop
    // - continues until allWithinTolerance returns true (specific number of individuals are within threshold)
    do {        
        // Genetic Crossover and mutation occur here
        //takes in oldAdults (the potential parents) and fills newAdults with descendants of the old adults
        newGeneration(oldAdults, newAdults, currentAnneal, generation, rng, cConstants);

        //fill oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected (numNans is for all adults in the generation - the oldAdults and the newAdults)
        preparePotentialParents(allAdults, newAdults, oldAdults, numNans, cConstants);

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

        //Perform utitlity tasks (adjusting anneal and reporting data)
        changeAnneal (newAdults, cConstants, new_anneal, currentAnneal, anneal_min, previousBestPosDiff, generation, posTolerance, dRate);
        reportGeneration (newAdults, cConstants, new_anneal, anneal_min, generation, numNans);

        // Before replacing new adults, determine whether all are within tolerance
        // Determines when loop is finished
        convergence = allWithinTolerance(posTolerance, speedTolerance, newAdults, cConstants);
        
        //Increment the generation counter
        ++generation;
    
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