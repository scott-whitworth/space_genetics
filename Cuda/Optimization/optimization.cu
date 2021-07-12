// Didymos Optimization Project using CUDA and a genetic algorithm

//TODO: Clarify complexities of the include paths
//TODO: What / why we are including
#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/individuals.h" // For individual structs, paths to rkParameters for randomParameters()
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // for testing rk4simple
#include "../Genetic_Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration()

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>
// Used to see if the best individual is changing when compared to a previous individual across generations 
// Returns true if the currentBest is not equal to previousBest within a distinguishable difference
// Input: previousBestPos - Position diference at the end of RK simulation, in AU
//        previousBestVel - Velocity difference, in AU/s (currently not implemented)
//        currentBest     - 'best' individual from current run, based on how individuals are sorted
//        distinguishRate - magnitude of difference
// Called inside of optimize to see if anneal rate needs to change

//! Currently not in use, changeInBest determined by cost, not posDiff
// bool changeInBest(double previousBestPos, double previousBestVel, const Individual & currentBest, double distinguishRate) {
//     //truncate is used here to compare doubles via the distinguguishRate, to ensure that there has been relatively no change.
//     if (trunc(previousBestPos/distinguishRate) != trunc(currentBest.posDiff/distinguishRate)) {
//         return true;
//     }
//     else {
//         /* //Used if Velocity should be considered
//         if (trunc(previousBestVel/distinguishRate) != trunc(currentBest.speedDiff/distinguishRate)) {
//             return true;
//         }
//         else return false;
//         */
//         return false;
//     }
// }

//Used to give rankings for sorting based on non-dominated sorting method.
//Assigns suitability rank to all individuals.
//MUST be called after cost has been assigned to all individuals (calling callRK)
//Input: pool - this generation of individuals, defined/initilized in optimimize
//       cConstants
void giveRank(Individual * pool, const cudaConstants* cConstants) {
    //non-denominated sorting method attempt
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

    // //Used to store the current front of individuals. first filled with the first front individuals(best out of all population)
    // std::vector<Individual> front;
    
    // //loop through each individual
    // for (int i = 0; i < cConstants->num_individuals; i++){
        
    //     //number of times pool[i] has been dominated
    //     pool[i].dominatedCount = 0;

    //     //set of solutions that pool[i] dominates
    //     // pool[i].dominated.clear();
    //     std::vector<Individual>().swap(pool[i].dominated);

    //     for(int j = 0; j < cConstants->num_individuals; j++){
            
    //         //if i dominates j, put j in the set of individuals dominated by i.
    //         if (dominates(pool[i], pool[j])){
    //             pool[i].dominated.push_back(pool[j]);
    //         }
    //         //if j dominates i, increase the number of times that i has been dominated
    //         else if (dominates(pool[j], pool[i])) {
    //             pool[i].dominatedCount++;
    //         }
    //     }
        
    //     //if i was never dominated, add it to the best front, front1. Making its ranking = 1.
    //     if (pool[i].dominatedCount == 0){
    //         pool[i].rank = 1;
    //         front.push_back(pool[i]);
    //     }
    //     std::cout << i << " ";
    // }
    // std::cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ";
    // //Used to assign rank number
    // int rankNum = 1;
    // //vector to store individuals in next front
    // std::vector<Individual> newFront;

    // //go until all individuals have been put in better ranks and none are left to give a ranking
    // while(!front.empty()) {
    //     //empty the new front to put new individuals in
    //     //newFront.clear();
    //     std::vector<Individual>().swap(newFront);

    //     //loop through all individuals in old front
    //     for(int i = 0; i < front.size(); i++){

    //         //loop through all the individuals that individual i dominated
    //         for(int j = 0; j < front[i].dominated.size(); j++){

    //             //subtract 1 from the dominated individuals' dominatedCount.
    //             //if an individual was dominated only once for example, it would be on the second front of individuals.
    //             front[i].dominated[j].dominatedCount--;

    //             //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
    //             if (front[i].dominated[j].dominatedCount == 0){
    //                 front[i].dominated[j].rank = i + 1;
    //                 newFront.push_back(front[i].dominated[j]);                        
    //             }
    //         }
    //     }
    //     //increment the rank number
    //     rankNum++;

    //     std::vector<Individual>().swap(front);
    //     //go to next front
    //     front = newFront;
    // }

    //Used to store the current front of individuals. first filled with the first front individuals(best out of all population)
    // filled with index of individuals in pool
    std::vector<int> front;
    
    //loop through each individual
    for (int i = 0; i < cConstants->num_individuals; i++){
        
        //number of times pool[i] has been dominated
        pool[i].dominatedCount = 0;

        //set of solutions that pool[i] dominates. Need to empty for each generation
        std::vector<int>().swap(pool[i].dominated);

        for(int j = 0; j < cConstants->num_individuals; j++){
            
            //if i dominates j, put the j index in the set of individuals dominated by i.
            if (dominates(pool[i], pool[j])){
                pool[i].dominated.push_back(j);
            }
            //if j dominates i, increase the number of times that i has been dominated
            else if (dominates(pool[j], pool[i])) {
                pool[i].dominatedCount++;
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (pool[i].dominatedCount == 0){
            pool[i].rank = 1;
            front.push_back(i);
        }
        //std::cout << i << " ";
    }
    //std::cout << "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ";
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

            //loop through all the individuals that the individual at index k dominated
            for(int l = 0; l < pool[front[k]].dominated.size(); l++){

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                pool[pool[front[k]].dominated[l]].dominatedCount--;

                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (pool[pool[front[k]].dominated[l]].dominatedCount == 0){

                    pool[pool[front[k]].dominated[l]].rank = k + 1;
                    
                    newFront.push_back(pool[front[k]].dominated[l]);                        
                }
            }
        }
        //increment the rank number
        rankNum++;

        std::vector<int>().swap(front);
        //go to next front
        front = newFront;
    }
    //for each individual, check if it is a NaN. If it is, give it a very low rank.
    for (int i = 0; i < cConstants->num_individuals; i++){
        if (pool[i].posDiff == 1.0){
            pool[i].rank = cConstants->num_individuals;
        }
    }
}

void giveDistance(Individual * pool, const cudaConstants* cConstants){

    for (int i = 0; i < cConstants->num_individuals; i++ ){
        //reset each individual's distance
        pool[i].distance = 0;
    }

    std::sort(pool, pool + cConstants->num_individuals, LowerPosDiff);
    pool[0].distance = 1.0e+12;
    pool[cConstants->num_individuals - 1].distance = 1.0e+12;

    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent solutions.
    for(int i = 1; i < cConstants->num_individuals - 1; i++){
        
        pool[i].distance = pool[i].distance + (pool[i+1].posDiff - pool[i-1].posDiff)/(pool[cConstants->num_individuals - 1].posDiff - pool[0].posDiff);
    }

    //Repeat above process for speedDiff    
    std::sort(pool, pool + cConstants->num_individuals, LowerSpeedDiff);
    pool[0].distance = 1.0e+12;
    pool[cConstants->num_individuals - 1].distance = 1.0e+12;

    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent solutions.
    for(int i = 1; i < cConstants->num_individuals - 1; i++){
        
        pool[i].distance = pool[i].distance + (pool[i+1].speedDiff - pool[i-1].speedDiff)/(pool[cConstants->num_individuals - 1].speedDiff - pool[0].speedDiff);
    }

}

bool changeInBest(double previousBestCost, const Individual & currentBest, double distinguishRate) {
    //truncate is used here to compare doubles via the distinguguishRate, to ensure that there has been relatively no change.
        if (trunc(previousBestCost/distinguishRate) != trunc(currentBest.cost/distinguishRate)) {
            return true;
        }
        else { 
            return false;
        }
}

// ** Assumes pool is sorted array of Individuals **
// Used in determining if main optimize loop continues
// Input: tolerance - posDiff threshold, determines max target distance
//        pool - this generation of Individuals, defined/initilized in optimimize
//        cConstants - struct holding config values, used for accessing best_count value
// Output: Returns true if top best_count individuals within the pool are within the tolerance
bool allWithinTolerance(double tolerance, Individual * pool, const cudaConstants* cConstants) {

    //Check what type of mission is running to use the correct cost function
    if (cConstants->missionType == Rendezvous){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            // if(pool[i].getCost_Soft(cConstants) >= tolerance) {
            //     //One was not within tolerance
            //     return false;
            // }
            if (pool[i].posDiff >= tolerance){
                return false;
            }

            if (pool[i].speedDiff >= tolerance){
                return false;
            }
        }
    }
    else if(cConstants->missionType == Impact){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            
            if(pool[i].getCost_Hard(cConstants) >= tolerance) {
                //One was not within 
                return false;
            }
        }  
    }

    // If iterated through and all were within tolerance, success
    return true;
}

// Main processing function for Genetic Algorithm
// - manages memory needs for genetic algorithm
// - deals with processing calls to CUDA callRK
// - exits when individuals converge on tolerance defined in Constants
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
    double absTol = cConstants->rk_tol; 
    // the starting step size for RK run
    // - note that the current step size varies throughout each run
    //TODO: Should this be based on max_numsteps?
    double stepSize = ((cConstants->orbitalPeriod) - timeInitial) / cConstants->GuessMaxPossibleSteps; 

    // Initial genetic anneal scalar
    double currentAnneal = cConstants->anneal_initial;

    // Main set of parameters for Genetic Algorithm
    // contains all thread unique input parameters
    Individual *inputParameters = new Individual[cConstants->num_individuals]; 

    // set to zero to force difference in first generation
    // double previousBestPos = 0; 
    // double previousBestVel = 0;
    double previousBestCost = 0;

    // Initilize individuals randomly or from a file
    if (cConstants->random_start) {
        // individuals set to randomly generated, but reasonable, parameters
        for (int i = 0; i < cConstants->num_individuals; i++) { 
            inputParameters[i] = Individual(randomParameters(rng, cConstants), cConstants);
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

            inputParameters[i] = Individual(example, cConstants);
        }
    }

    // Collection of individuals used in the genetic selection process
    //  - filled in selectSurvivors, based on callRK output
    //  - stores the winners of the head-to-head competition
    Individual *survivors = new Individual[cConstants->survivor_count]; 

    // Number of individuals that need to be evaluated
    // - the whole population is in first loop
    // - subsequent generations only calculate *new* individuals
    int newInd = cConstants->num_individuals;

    // number of current generation
    double generation = 0;    
    
    // how far away the best individual is from the tolerance value
    double currentCost; 
    // Genetic solution tolerance 
    // - (currently just the position threshold which is furthest distance from the target allowed)
    // - could eventually take into account velocity too and become a more complex calculation
    double tolerance = cConstants->pos_threshold; 
             
    // distinguishable rate used in changeInBest()
    //  - used to help check for a change in anneal
    //  - Gets smaller when no change is detected
    double dRate = 1.0e-8;

    // Flag for finishing the genetic process
    // set by allWithinTolerance()
    bool convergence = false;

    // main gentic algorithm loop
    // - continues until allWithinTolerance returns true (specific number of individuals are within threshold)
    do {
        // each inputParameter represents an individual set of starting parameters
        // GPU based runge kutta process determines final position and velocity based on parameters
        // newInd - how many individuals that are *new* that need to be evaluated
        //        - All individuals first generation
        //        - only new individuals, from crossover, in subsequent generations
        // (inputParameters + (cConstants->num_individuals - newInd)) value accesses the start of the section of the inputParameters array that contains new individuals
        callRK(newInd, cConstants->thread_block_size, inputParameters + (cConstants->num_individuals - newInd), timeInitial, stepSize, absTol, calcPerS, cConstants); // calculate trajectories for new individuals

        //numNans - number of times a nan is found in 100 generations.
        int numNans = 0;

        // if we got bad results reset the Individual to random starting values (it may still be used for crossover) and set the final position to be way off so it gets replaced by a new Individual
        for (int k = 0; k < cConstants->num_individuals; k++) {
            //Checking each individuals final position for NaNs
            if (isnan(inputParameters[k].finalPos.r) || isnan(inputParameters[k].finalPos.theta) || isnan(inputParameters[k].finalPos.z) || isnan(inputParameters[k].finalPos.vr) || isnan(inputParameters[k].finalPos.vtheta) || isnan(inputParameters[k].finalPos.vz)) {
                //std::cout << std::endl << std::endl << "NAN FOUND" << std::endl << std::endl;
                numNans++;
                inputParameters[k] = Individual(randomParameters(rng, cConstants), cConstants);
                // Set to be a bad individual by giving it bad posDiff and speedDiffs
                // therefore also having a bad cost value
                // won't be promoted in crossover
                inputParameters[k].posDiff = 1.0;

                if (cConstants->missionType == Rendezvous){
                    inputParameters[k].speedDiff = 1.0;//This is an undesirable result for an rendezvous mission (approx. 500c!)
                    // calculate its new cost function based on 'bad' differences
                    inputParameters[k].getCost_Soft(cConstants);
                }
                else if (cConstants->missionType == Impact){
                    inputParameters[k].speedDiff = 0.0;//This is an undesirable result for an impact mission
                    // calculate its new cost function based on 'bad' differences
                    inputParameters[k].getCost_Hard(cConstants);                   
                }
            }
            
            //if (k > 0 && (inputParameters[k].posDiff != 1.0)){
            //if (inputParameters[k].posDiff != 1.0){
                //calculate the % difference of the individual to the best cost individual
                //inputParameters[k].difference = checkDifference(inputParameters[0],inputParameters[k]); 

                //if(inputParameters[0].posDiff > 1.0e-06){

                    // if(inputParameters[k].difference < 20.0){
                    //     inputParameters[k].cost = inputParameters[k].cost * 100000;
                    // }
                    // else if (inputParameters[k].difference < 50.0 && inputParameters[k].difference >= 20.0){
                    //     inputParameters[k].cost = inputParameters[k].cost * 10000;
                    // }
                    // else if (inputParameters[k].difference < 70.0 && inputParameters[k].difference >= 50.0){
                    //     inputParameters[k].cost = inputParameters[k].cost * 1000;
                    // }
                    // else if (inputParameters[k].difference < 90.0 && inputParameters[k].difference >= 70.0){
                    //     inputParameters[k].cost = inputParameters[k].cost * 100;
                    // }

                    //low speed penalty
                    // if(inputParameters[k].speedDiff > 1.0E-09){
                    //     inputParameters[k].cost = inputParameters[k].cost * 10000;
                    // }
                    // else if (inputParameters[k].speedDiff < 1.0E-08 && inputParameters[k].speedDiff >= 1.0E-09){
                    //     inputParameters[k].cost = inputParameters[k].cost * 1000;
                    // }

                //}
                // else {
                //     if(inputParameters[k].difference < 80.0){
                //         inputParameters[k].cost = inputParameters[k].cost * 100000;
                //     }
                //     // else if (inputParameters[k].difference < 0.1 && inputParameters[k].difference >= 0.01){
                //     //     inputParameters[k].cost = inputParameters[k].cost * 10000;
                //     // }
                //     // else if (inputParameters[k].difference < 1.0 && inputParameters[k].difference >= 0.1){
                //     //     inputParameters[k].cost = inputParameters[k].cost * 1000;
                //     // }
                //     // else if (inputParameters[k].difference < 10.0 && inputParameters[k].difference >= 1.0){
                //     //     inputParameters[k].cost = inputParameters[k].cost * 100;
                //     // }
                // }
            //}
        }

        // Preparing survivor pool with individuals for the newGeneration crossover
        // Survivor pool contains:
        //               - individuals with best PosDiff
        //               - individuals with best speedDiffs
        //               - depends on cConstants->sortingRatio (0.1 is 10% are best PosDiff for example)
        // inputParameters is left sorted by individuals with best speedDiffs 
        selectSurvivors(inputParameters, cConstants->num_individuals, cConstants->survivor_count, survivors, cConstants->sortingRatio, cConstants->missionType); // Choose which individuals are in survivors, current method selects half to be best posDiff and other half to be best speedDiff

        // sort individuals based on overloaded relational operators
        // gives reference of which to replace and which to carry to the next generation
        giveRank(inputParameters, cConstants);
        giveDistance(inputParameters, cConstants);
        std::sort(inputParameters, inputParameters + cConstants->num_individuals, rankDistanceSort);

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the terminalDisplay() calls across generations 
        std::cout << '.';

        // Calculate how far best individual is from the ideal cost value (currently is the positionalDifference of the best individual)
        // TODO: Change this later to take into account more than just the best individual and its position difference
        currentCost = inputParameters[0].cost; 

        // Scaling anneal based on proximity to tolerance
        // Far away: larger anneal scale, close: smaller anneal
        double new_anneal = currentAnneal * (1 - tolerance / currentCost);
        
        //Process to see if anneal needs to be adjusted
        // If generations are stale, anneal drops
        Individual currentBest;
        // Compare current best individual to that from CHANGE_CHECK many generations ago.
        // If they are the same, change size of mutations
        if (static_cast<int>(generation) % cConstants->change_check == 0) { 
            currentBest = inputParameters[0];
            // checks for anneal to change
            // previousBest starts at 0 to ensure changeInBest = true on generation 0
            if ( !(changeInBest(previousBestCost, currentBest, dRate)) ) { 
                //this ensures that changeInBest never compares two zeros, thus keeping dRate in relevance as the posDiff lowers
                if (trunc(currentBest.cost/dRate) == 0) { 
                    while (trunc(currentBest.cost/dRate) == 0) {
                        dRate = dRate/10; 
                    }
                    std::cout << "\nnew dRate: " << dRate << std::endl;
                }
                // If no change in BestIndividual across generations, multiply currentAnneal with anneal factor
                currentAnneal = currentAnneal * cConstants->anneal_factor;
                std::cout << "\nnew anneal: " << currentAnneal << std::endl;
            }

            // previousBestPos = currentBest.posDiff;
            // previousBestVel = currentBest.speedDiff;
            previousBestCost = currentBest.cost;
        }

        // If in recording mode and write_freq reached, call the record method
        if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
            std::sort(inputParameters, inputParameters + cConstants->num_individuals);
            recordGenerationPerformance(cConstants, inputParameters, generation, new_anneal, cConstants->num_individuals);
        }
        std::sort(inputParameters, inputParameters + cConstants->num_individuals, rankDistanceSort);
        // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
            // Prints the best individual's posDiff / speedDiff and cost
            //terminalDisplay(inputParameters[0], generation);

            //best position individual
            std::cout << "\nBest Position Individual: \n";
            std::sort(inputParameters, inputParameters + cConstants->num_individuals, LowerPosDiff);
            terminalDisplay(inputParameters[0], generation);

            if(cConstants->missionType == Rendezvous){
                //Best lower speed individual
                std::cout << "\nBest Speed Individual: \n";
                std::sort(inputParameters, inputParameters + cConstants->num_individuals, LowerSpeedDiff);
                terminalDisplay(inputParameters[0], generation);
            }
            else if(cConstants->missionType == Impact){
                //Best higher speed individual
                std::cout << "\nBest Speed Individual: \n";
                std::sort(inputParameters, inputParameters + cConstants->num_individuals, HigherSpeedDiff);
                terminalDisplay(inputParameters[0], generation);
            }

            //reset array
            std::sort(inputParameters, inputParameters + cConstants->num_individuals);
            terminalDisplay(inputParameters[0], generation);
            std::cout << "\n# of Nans this increment: " << numNans << "\n" << std::endl;
            numNans = 0; //Reset the tally of nans.

            recordAllIndividuals(cConstants, inputParameters, cConstants->num_individuals, generation);
        }
        std::sort(inputParameters, inputParameters + cConstants->num_individuals, rankDistanceSort);
        // Before replacing new individuals, determine whether all are within tolerance
        // Determines when loop is finished
        convergence = allWithinTolerance(tolerance, inputParameters, cConstants);

        // Create a new generation and increment the generation counter
        // Genetic Crossover and mutation occur here
        newInd = newGeneration(survivors, inputParameters, cConstants->survivor_count, cConstants->num_individuals, new_anneal, cConstants, rng, generation);
        ++generation;
    
        //Loop exits based on result of allWithinTolerance and if max_generations has been hit
    } while ( !convergence && generation < cConstants->max_generations);

    // Call record for final generation regardless of frequency
    // for the annealing argument, set to -1 (since the anneal is only relevant to the next generation and so means nothing for the last one)
    if (cConstants->record_mode == true) {
        recordGenerationPerformance(cConstants, inputParameters, generation, -1, cConstants->num_individuals);
    }
    // Only call finalRecord if the results actually converged on a solution
    // also display last generation onto terminal
    if (convergence) {
        terminalDisplay(inputParameters[0], generation);
        finalRecord(cConstants, inputParameters, static_cast<int>(generation));
    }
    
    delete [] inputParameters;
    delete [] survivors;

    return calcPerS;
}

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
