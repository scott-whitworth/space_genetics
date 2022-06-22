// Didymos Optimization Project using CUDA and a genetic algorithm

//TODO: Clarify complexities of the include paths
//TODO: What / why we are including
#include "../Earth_calculations/earthInfo.h"  // For launchCon and EarthInfo()
#include "../Genetic_Algorithm/adult.h" // For adult structs, paths to rkParameters for randomParameters()
#include "../Genetic_Algorithm/child.h" // For child structs, paths to rkParameters for randomParameters()
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // for testing rk4simple
#include "../Genetic_Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration()
#include "../Genetic_Algorithm/genetic_algorithm.h" //For functions that set up new generations
#include "../Genetic_Algorithm/sort.h" //For functions that will allow for sorting of the adult arrays by giving them ranks and distances
#include "../Genetic_Algorithm/anneal.h" //For all the annealing functions
//
//#include "../Unit_Testing/testing_sorts.cpp"

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>   // allows us to use vectors instead of just arrays

//----------------------------------------------------------------------------------------------------------------------------
// ** Assumes pool is sorted array of Adults **
// Used in determining if main optimize loop continues
// Input: posTolerance - posDiff threshold, determines max target distance
//        speedTolerance - speedDiff threshold, determines max speed at target
//        oldAdults - this generation of Adults, defined/initialized in optimimize
//        cConstants - struct holding config values, used for accessing best_count value
// Output: Returns true if top best_count adults within the pool are within the tolerance
bool checkTolerance(double posTolerance, double speedTolerance, const std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// TEST / LIKELY TEMPORARY FUNCTION
// This function will find the minimum, maximum, and average distance, the avg age, and the avg and max birthdays of the individuals in allAdults, which will then be used for reporting
// 
// Inputs:  allAdults - array of adults that will be considered\
//          generation - the current generation
//          minDist - minimum distance that will be calculated
//          avgDist - the average distance that will be calculated
//          maxDist - the maximum distance that will be calculated
//          avgAge  - the avegrage age of the adults, relative to the current generation
//          avgBirthday - average birth generation for the adults
//          oldestBirthday - the oldest adult's birth generation
// Outputs: The arguments will be filled in with the up-to-date values for this generation
void calculateGenerationValues (const std::vector<Adult> & allAdults, const int & generation, double & minDist, double & avgDist, double & maxDist, double & avgAge, double & avgBirthday, int & oldestBirthday);

//----------------------------------------------------------------------------------------------------------------------------
//Input: all the updated parameters of the current generation
//Output: calls the various terminal display functions when needed
//Function that handles the reporting of a generation's performance
void reportGeneration (std::vector<Adult>& oldAdults, std::vector<Adult> & allAdults, const cudaConstants* cConstants, const double & currentAnneal, const double & anneal_min, const int & generation, int & numNans);

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

//Returns true if top best_count adults within the oldAdults vector are within the tolerance
bool checkTolerance(double posTolerance, double speedTolerance, const std::vector<Adult>& oldAdults, const cudaConstants* cConstants) {

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

//Function that will calculate distance and birthday values for a generation
void calculateGenerationValues (const std::vector<Adult> & allAdults, const int & generation, double & minDist, double & avgDist, double & maxDist, double & avgAge, double & avgBirthday, int & oldestBirthday){
    //Reset the dist values
    minDist = 2; //Set the min dist to the maximum possible value, so that it will be changed
    avgDist = 0; 
    maxDist = 0; //Set the max dist to the min possible value, so that it is garunteed to be changed

    //Reset the age values
    avgAge = 0;
    avgBirthday = 0; 
    oldestBirthday = generation; //Set to generation, since it is the "newest" possible oldest birthday
    
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

        //Check to see if this adult's birthday is older than the current oldest
        if (allAdults[i].birthday < oldestBirthday) {
            oldestBirthday = allAdults[i].birthday; 
        }

        //Add to the avg distance
        avgDist += allAdults[i].distance;

        //Add to the avg age values
        //avgAge measures the age relative to the current generation (how many generations old), so it is generation minus the adults' birthdays 
        avgBirthday += allAdults[i].birthday;
        avgAge += (generation - allAdults[i].birthday);   
    }

    //Divide the averages by the number of adults to get the averages
    avgDist /= allAdults.size();
    avgAge /= allAdults.size(); 
    avgBirthday /= allAdults.size(); 
    
}

//Function that writes the results of the inserted generation
void reportGeneration (std::vector<Adult> & oldAdults, std::vector<Adult> & allAdults, const cudaConstants* cConstants, const double & currentAnneal, const double & anneal_min, const int & generation, int & numNans){
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
        recordAllIndividuals("NextParents", cConstants, allAdults, allAdults.size(), generation);
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
    //double timeInitial = 0; // the starting time of the trip is always defined as zero   

    // Runge Kutta adaptive time step error tolerance
    // Not used now that callRK has been moved to ga_crossover
    //double absTol = cConstants->rk_tol; 

    // Initial genetic anneal scalar
    double currentAnneal = cConstants->anneal_initial;

    //lower bound for anneal so it does not get too small. Only used with rendezvous mission.
    double anneal_min = cConstants->anneal_initial;
    
    // Main set of parameters for Genetic Algorithm
    // contains all thread unique input parameters
    // The "children pool" of the current genertation
    //TODO: Where is this size set? Are we certian the size is correct? / do we have positive control?
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
    //TODO: Figure out if we can delete

    // number of current generation
    int generation = 0;    
    //TODO: This, for sure, should be an int
    //      This will be a major refactor, do it very carefully.
    
    // Genetic solution tolerance for each objective
    double posTolerance = cConstants->pos_threshold;
    double speedTolerance = cConstants->speed_threshold;  
    //TODO: Can we just use the cConstants value directly? I don't think we need the doubles defined here
             
    // distinguishable rate used in changeInBest()
    //  - used to help check for a change in anneal
    //  - Gets smaller when no change is detected
    double dRate = 1.0e-8;

    // Flag for finishing the genetic process
    // set by checkTolerance()
    bool convergence = false;

    //sets the anneal for the generation
    //used for proximity based annealing
    double new_anneal;

    //number of Nans, specifically used for diagnostic recording
    //  couting all adults in the generation - includes oldAdults and newAdults that have nan values
    int numNans = 0;

    //Initialize variables needed for distance and birthday reporting
    double maxDistance, minDistance, avgDistance, avgAge, avgBirthday;
    int oldestBirthday;

    //Creates the individuals needed for the 0th generation
    //Need to make children, then callRK, then make into adults (not currently doing that)
    //oldAdults goes into this function empty and comes out filled with num_individuals Adults
    //      these adults are either randomly generated or pulled from a file
    createFirstGeneration(oldAdults, cConstants, rng, generation); 

    // main gentic algorithm loop
    // - continues until checkTolerance returns true (specific number of individuals are within threshold)
    do {        
        // Genetic Crossover and mutation occur here
        //takes in oldAdults (the potential parents) and fills newAdults with descendants of the old adults
        //oldAdults is filled with the potential parents for a generation (num_individuals size) 
        //      after the run, oldAdults should remain the same
        //newAdults is empty and will be filled with the "grown" children generated in this method
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE NEW GEN-_-_-_-_-_-_-_-_-_\n\n";
        newGeneration(oldAdults, newAdults, currentAnneal, generation, rng, cConstants);
        //TODO: What is the state of oldAdults / newAdults after this is run (just a reminder)

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE PREP PARENTS-_-_-_-_-_-_-_-_-_\n\n";
        //fill oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected (numNans is for all adults in the generation - the oldAdults and the newAdults)
        //allAdults will be filled with the last generation's set of parents and their offspring when this starts (sorted by rankDistanceSort)
        //      by the end of this function, it will be filled with the new generation's set of parents and children sorted by rankDistanceSort
        //newAdults goes in with the "grown" children created in new generation (size num_individuals)
        //      by the end of the function, it is cleared
        //oldAdults goes in with the pool of potential parents that may have generated the newAdults
        //      by the end of the function, it is filled with the best num_individuals adults from allAdults (sorted by rankDistanceSort) 
        preparePotentialParents(allAdults, newAdults, oldAdults, numNans, cConstants, generation);
        //TODO: What is the state of all/old/newAdults (just a reminder)


        //UPDATE (6/7/22): callRK moved into ga_crossover, the effect on the efficiency of the code is not known yet
        //Output survivors for the current generation if write_freq is reached
        if (generation % cConstants->all_write_freq == 0 && cConstants->record_mode == true) {
            recordAllIndividuals("Survivors", cConstants, oldAdults, cConstants->survivor_count, generation);
        }

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the terminalDisplay() calls across generations 
        std::cout << '.';

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE ANNEAL STFF-_-_-_-_-_-_-_-_-_\n\n";
        //Perform utitlity tasks (adjusting anneal and reporting data)
        //Assumes oldAdults is in rankDistance order
        changeAnneal (oldAdults, cConstants, new_anneal, currentAnneal, anneal_min, previousBestPosDiff, generation, posTolerance, dRate);

        //Calculate variables for birthdays and distances
        calculateGenerationValues(allAdults, generation, minDistance, avgDistance, maxDistance, avgAge, avgBirthday, oldestBirthday);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE RECORD-_-_-_-_-_-_-_-_-_\n\n";
        reportGeneration (oldAdults, allAdults, cConstants, currentAnneal, anneal_min, generation, numNans);
        //TODO: verify that old/AllAdults are in the right order (and clarify with comment here)

        // Before replacing new adults, determine whether all are within tolerance
        // Determines when loop is finished
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CONVERGENCE CHECK-_-_-_-_-_-_-_-_-_\n\n";
        convergence = checkTolerance(posTolerance, speedTolerance, oldAdults, cConstants);
        
        //Increment the generation counter
        ++generation;

        //sorts oldAdults using rankDistance sort
        //This will prep the oldAdults array to be used to generate new children; this makes sure any sorts from the record/report functions are undone
        std::sort(oldAdults.begin(), oldAdults.end(), rankDistanceSort); 
        //TODO: If we need this for generateNewChildren, why not just put this sort in generateNewChildren?
    
        //Loop exits based on result of checkTolerance and if max_generations has been hit
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
        finalRecord(cConstants, oldAdults, generation);
    }

    return calcPerS;
}