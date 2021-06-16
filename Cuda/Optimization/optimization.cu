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
            if(pool[i].getCost_Soft(cConstants) >= tolerance) {
                //One was not within tolerance
                return false;
            }
        }
    }
    else if(cConstants->missionType == Impact){
        // Iterate to check best_count number of 'top' individuals
        for (int i = 0; i < cConstants->best_count; i++) {
            if(pool[i].getCost_Hard(cConstants) >= tolerance) {
                //One was not within tolerance
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
<<<<<<< HEAD
                inputParameters[k].speedDiff = 1.0;
=======
>>>>>>> 88fdfdef3bbc4974891fbac0193337bbfa1180a6

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
        std::sort(inputParameters, inputParameters + cConstants->num_individuals);

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
                if (trunc(currentBest.posDiff/dRate) == 0) { 
                    while (trunc(currentBest.posDiff/dRate) == 0) {
                        dRate = dRate/10; 
                    }
                    std::cout << "\nnew dRate: " << dRate << std::endl;
                }
                // If no change, multiply currentAnneal with anneal factor
                currentAnneal = currentAnneal * cConstants->anneal_factor;
                std::cout << "\nnew anneal: " << currentAnneal << std::endl;
            }

            // previousBestPos = currentBest.posDiff;
            // previousBestVel = currentBest.speedDiff;
            previousBestCost = currentBest.cost;
        }

        // If in recording mode and write_freq reached, call the record method
        if (static_cast<int>(generation) % cConstants->write_freq == 0 && cConstants->record_mode == true) {
            recordGenerationPerformance(cConstants, inputParameters, generation, new_anneal, cConstants->num_individuals);
        }
        
        // Only call terminalDisplay every DISP_FREQ, not every single generation
        if ( static_cast<int>(generation) % cConstants->disp_freq == 0) {
            // Prints the best individual's posDiff / speedDiff and cost
            //terminalDisplay(inputParameters[0], generation);

            //best position individual
            std::cout << "\nBest Position Individual: \n";
            std::sort(inputParameters, inputParameters + cConstants->num_individuals, LowerPosDiff);
            terminalDisplay(inputParameters[0], generation);

            //Best speed individual
            std::cout << "\nBest Speed Individual: \n";
            std::sort(inputParameters, inputParameters + cConstants->num_individuals, LowerSpeedDiff);
            terminalDisplay(inputParameters[0], generation);

            //reset array
            std::sort(inputParameters, inputParameters + cConstants->num_individuals);
            terminalDisplay(inputParameters[0], generation);
            std::cout << "\n# of Nans this increment: " << numNans << "\n" << std::endl;
            numNans = 0; //Reset the tally of nans.
        }

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
