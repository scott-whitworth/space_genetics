// Spaceflight Optimization Project using CUDA and a genetic algorithm

#include "../Planet_calculations/planetInfo.h"  // For launchCon and EarthInfo(); includes elements, runge_kutta, & config.h
#include "../Genetic_Algorithm/adult.h" // For adult structs, paths to rkParameters for randomParameters(); includes child, rkParameters, & planetInfo
#include "../Genetic_Algorithm/child.h" // For child structs, paths to rkParameters for randomParameters(); includes config.h, rkParameters, & planetInfo
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" // for testing rk4simple; includes calcFourier, motion_eqns, child, & gpuMem
#include "../Genetic_Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration(); includes constants.h, adult, & child
#include "../Genetic_Algorithm/genetic_algorithm.h" //For functions that set up new generations; includes constants.h, adult, child, ga_crossover, & sort
#include "../Genetic_Algorithm/referencePoints.h" //For the ReferencePoints class which deals with calculating reference points and setting adult's rarity
#include "../Genetic_Algorithm/sort.h" //For functions that will allow for sorting of the adult arrays by giving them ranks and distances; includes constants.h & adult
#include "../Genetic_Algorithm/anneal.h" //For all the annealing functions; includes constants.h & adult
#include "../Runge_Kutta/gpuMem.cuh" // for initializing and deallocating; includes child, rkParameters, & config.h

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>   // allows us to use vectors instead of just arrays
#include <string>
#include <chrono>
#include <algorithm>


//----------------------------------------------------------------------------------------------------------------------------
// ** Assumes pool is sorted array of Adults **
// Used in determining if main optimize loop continues
// Input: oldAdults - this generation of Adults, defined/initialized in optimimize
//        cConstants - struct holding config values, used for accessing best_count value and objectives list
// Output: Returns true if top best_count adults within the pool are within the tolerance
bool checkTolerance(std::vector<Adult> & oldAdults, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// TEST / LIKELY TEMPORARY FUNCTION
// This function will find the minimum, maximum, and average distance, average pos and speed diffs, the number of duplicates,
//                         the avg age, and the avg and max birthdays of the individuals in allAdults, which will then be used for reporting
// 
// Inputs:  allAdults - array of adults that will be considered
//          objectives - the vector of this run's objectives
//          objectiveAvgValues - a vector which will hold the generations's average parameter values for each of the objectives
//          duplicateNum - the number of duplicate adults found
//          minDist - minimum distance that will be calculated
//          avgDist - the average distance that will be calculated
//          maxDist - the maximum distance that will be calculated
//          generation - the current generation
//          avgAge  - the avegrage age of the adults, relative to the current generation
//          avgBirthday - average birth generation for the adults
//          oldestBirthday - the oldest adult's birth generation
// Outputs: The arguments will be filled in with the up-to-date values for this generation
void calculateGenerationValues (const std::vector<Adult> & allAdults, const std::vector<objective> & objectives, std::vector<double> & objectiveAvgValues, int & duplicateNum, double & minDist, double & avgDist, double & maxDist, int & minSteps, int & avgSteps, int & maxSteps, const int & generation, double & avgAge, double & avgBirthday, int & oldestBirthday);

//----------------------------------------------------------------------------------------------------------------------------
// Main processing function for Genetic Algorithm
// - manages memory needs for genetic algorithm
// - deals with processing calls to CUDA callRK
// - exits when individuals converge on tolerance defined in Constants
double optimize(const cudaConstants* cConstants, GPUMem & gpuValues, const ReferencePoints & refPoints);

//-----------------------------------------------------------------------------------------------------------------------------
int main () {
    
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "\n\nDevice Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    cudaSetDevice(0);

    // Declare the genetic constants used, with file path being used to receive initial values
    cudaConstants * cConstants = new cudaConstants(); 

    //preallocates all the memory for the varaibles used by the GPU
    //also allows the GPU to access the marsLaunchCon without reloading it everytime
    GPUMem gpuValues;

    //Creates the reference points for the rest of the program
    ReferencePoints refPoints(cConstants);

    //Display the number of reference points created as a sanity check
    std::cout << "\n" << refPoints.points.size() << " reference points created.\n";

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

        // pre-calculate a table of Earth's and Mars' position within possible mission time range
        // defined as global variable
        // accessed on the CPU when individuals are initialized
        launchCon = new PlanetInfo(cConstants, EARTH); 
        marsLaunchCon = new PlanetInfo(cConstants, MARS);
        //This ensures that we copy the correct size of marsCon to the GPU
        int marsConSize = getPlanetSize(cConstants);
        //initialize all values needed for GPU calculations
        gpuValues.initialize(cConstants, marsConSize, marsLaunchCon->getAllPositions());

        // Call optimize with the current parameters in cConstants
        optimize(cConstants, gpuValues, refPoints);
        
        // Deallocate launchCon info for this run as it may be using a different time range in the next run
        delete launchCon; 
        delete marsLaunchCon;
        gpuValues.free();
    }
    // Now that the optimize function is done (assumed that optimize() also records it), deallocate memory of the cudaConstants
    delete cConstants;
    
    return 0;
}

//Returns true if top best_count adults within the oldAdults vector are within the tolerance
bool checkTolerance(std::vector<Adult>& oldAdults, const cudaConstants* cConstants) {
    //Sort the vector by rank distance to make sure the program checks the correct adult
    std::sort(oldAdults.begin(), oldAdults.end(), rankRaritySort); 

    //The function needs to check if the best adult meets the convergence tolerance for each objective
    //Iterate through the objectives
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //Check to see if the top best_count adults have met convergence for this parameter
        for (int j = 0; j < cConstants->best_count; j++) {

            //Check to see if the adult's parameter is larger than the convergence 
            if (oldAdults[j].getParameters(cConstants->missionObjectives[i]) > cConstants->missionObjectives[i].convergenceThreshold) {
                //Return false as a parameter that needs to be minimized is larger than the convergence threshold
                return false;
            }
        }

        // if (cConstants->missionObjectives[i].goal < 0) {//Minimization
        //     //Check to see if the top best_count adults have met convergence for this parameter
        //     for (int j = 0; j < cConstants->best_count; j++) {

        //         //Check to see if the adult's parameter is larger than the convergence 
        //         if (oldAdults[j].getParameters(cConstants->missionObjectives[i]) > cConstants->missionObjectives[i].convergenceThreshold) {
        //             //Return false as a parameter that needs to be minimized is larger than the convergence threshold
        //             return false;
        //         }
        //     }
        // }
        // else if (cConstants->missionObjectives[i].goal > 0) {//Maximization
        //     //Check to see if the top best_count adults have met convergence for this parameter
        //     for (int j = 0; j < cConstants->best_count; j++) {
        //         //Check to see if the adult's parameter is smaller than the convergence 
        //         if (oldAdults[j].getParameters(cConstants->missionObjectives[i]) < cConstants->missionObjectives[i].convergenceThreshold) {
        //             //Return false as a parameter that needs to be maximized is smaller than the convergence threshold
        //             return false;
        //         }
        //     }
        // }
        // //No mission type was identified 
        // else {
        //     std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        // }
    }

    //If the program reaches this spot, it means all of the adult's parameters have met the convergence threshold
    //  Otherwise, the function would have already returned false
    //  Thus, the adult has converged and it is appropriate to return true
    return true; 
}

//Function that will calculate distance and birthday values for a generation
void calculateGenerationValues (const std::vector<Adult> & allAdults, const std::vector<objective> & objectives, std::vector<double> & objectiveAvgValues, int & duplicateNum, double & minDist, double & avgDist, double & maxDist, int & minSteps, int & avgSteps, int & maxSteps, const int & generation, double & avgAge, double & avgBirthday, int & oldestBirthday){
    //Reset the dist values
    minDist = 2; //Set the min dist to the maximum possible value, so that it will be changed
    avgDist = 0; 
    maxDist = 0; //Set the max dist to the min possible value, so that it is garunteed to be changed

    minSteps = 1000000;
    avgSteps = 0;
    maxSteps = 0;

    //Reset the average parameter values
    objectiveAvgValues.clear();
    //Prep the avg value vector to have an index for each objective
    for (int i = 0; i < objectives.size(); i++) {
        objectiveAvgValues.push_back(0.0); 
    }
    
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
        else if (allAdults[i].distance > maxDist) {
            maxDist = allAdults[i].distance; //Ser the new max distance
        }

        //Check to see if this adult's step count is a new minimum
        if (allAdults[i].stepCount < minSteps) {
            minSteps = allAdults[i].stepCount; //Set the new min distance
        }
        //Check to see if this adult's step count is the new maximum
        else if (allAdults[i].stepCount > maxSteps) {
            maxSteps = allAdults[i].stepCount; //Set the new max distance
        }

        //Check to see if this adult's birthday is older than the current oldest
        if (allAdults[i].birthday < oldestBirthday) {
            oldestBirthday = allAdults[i].birthday; 
        }
           
        //Add to the avg distance
        avgDist += allAdults[i].distance;
        avgSteps += allAdults[i].stepCount;

        //Add the adult's parameter values to the necessary spot in the objective average value vector
        for (int j = 0; j < objectives.size(); j++) {
            //See if it is a minimization or maximization objective
            //If it is a maximization the individual's total needs to be subtracted from the avg value as the value is stored as a negative in the individual
            if (objectives[j].goal < 0) {
                //Minimization
                objectiveAvgValues[j] += allAdults[i].getParameters(objectives[j]);
            }
            else {
                //Maximization
                objectiveAvgValues[j] -= allAdults[i].getParameters(objectives[j]);
            }
            
        }
        
        //Add to the avg age values
        //avgAge measures the age relative to the current generation (how many generations old), so it is generation minus the adults' birthdays 
        avgBirthday += allAdults[i].birthday;
        avgAge += (generation - allAdults[i].birthday);   
    }
    //Possible floating point roundoff error
    //  Shouldn't matter, as we don't ever compared averages to individuals
    //Divide the averages by the number of adults to get the averages
    avgDist /= allAdults.size();
    avgSteps /= allAdults.size();
    for (int i = 0; i < objectiveAvgValues.size(); i++) {
        objectiveAvgValues[i] /= allAdults.size();
    }
    avgAge /= allAdults.size(); 
    avgBirthday /= allAdults.size(); 
    
}

//Function that will facilitate the process of finding an optimal flight path
double optimize(const cudaConstants* cConstants, GPUMem & gpuValues, const ReferencePoints & refPoints) {
    // Not used, previously used for reporting computational performance
    double calcPerS = 0;

    time_t timeSeed = cConstants->time_seed;
    std::mt19937_64 rng(timeSeed); // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;

    //Initialize the output object with the base folder location (..\Output_Files\)
    output genOutputs(cConstants);

    // Initial genetic anneal scalar
    double currentAnneal = cConstants->anneal_initial;
    
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

    // number of current generation
    int generation = 0;

    genOutputs.recordMarsData(cConstants,generation);
    genOutputs.recordEarthData(cConstants,generation);

    // Flag for finishing the genetic process
    // set by checkTolerance()
    bool convergence = false;

    //number of errors, specifically used for diagnostic recording
    //  couting all adults in the generation - includes oldAdults and newAdults that have nan values
    int numErrors = 0;
    int marsErrors = 0;

    //Initialize variables needed for distance, number of duplicate adults, and birthday reporting
    int duplicateNum = 0;
    double maxDistance, minDistance, avgDistance, avgAge, avgBirthday;
    int minSteps, avgSteps, maxSteps, oldestBirthday;

    //Inititalize variables for storing the average time per generation
    std::chrono::time_point<std::chrono::system_clock> runStartTime = std::chrono::system_clock::now();
    std::chrono::duration<float> totRunTime;

    //double worstOPD = 0;

    //Vector used to report the average parameter value for each objective
    std::vector<double> objectiveAvgValues; 

    //Creates the individuals needed for the 0th generation
    //Need to make children, then callRK, then make into adults (not currently doing that)
    //oldAdults goes into this function empty and comes out filled with num_individuals Adults
    //      these adults are either randomly generated or pulled from a file
    createFirstGeneration(oldAdults, cConstants, rng, generation, gpuValues, refPoints); 


    // main gentic algorithm loop
    // - continues until checkTolerance returns true (specific number of individuals are within threshold)
    do {
        // Genetic Crossover and mutation occur here
        //takes in oldAdults (the potential parents) and fills newAdults with descendants of the old adults
        //oldAdults is filled with the potential parents for a generation (num_individuals size) 
        //      after the run, oldAdults should remain the same
        //newAdults is empty and will be filled with the "grown" children generated in this method
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE NEW GEN-_-_-_-_-_-_-_-_-_\n\n";
        newGeneration(oldAdults, newAdults, currentAnneal, generation, rng, cConstants, gpuValues);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE PREP PARENTS-_-_-_-_-_-_-_-_-_\n\n";
        //fill oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected (numErrors is for all adults in the generation - the oldAdults and the newAdults)
        //allAdults will be filled with the last generation's set of parents and their offspring when this starts (sorted by rankDistanceSort)
        //      by the end of this function, it will be filled with the new generation's set of parents and children sorted by rankDistanceSort
        //newAdults goes in with the "grown" children created in new generation (size num_individuals)
        //      by the end of the function, it is cleared
        //oldAdults goes in with the pool of potential parents that may have generated the newAdults
        //      by the end of the function, it is filled with the best num_individuals adults from allAdults (sorted by rankDistanceSort) 
        preparePotentialParents(allAdults, newAdults, oldAdults, numErrors, duplicateNum, cConstants, refPoints, generation, currentAnneal, marsErrors);

        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the terminalDisplay() calls across generations 
        std::cout << '.';

        //Perform utitlity tasks (adjusting anneal and reporting data)
        //Calculate variables for birthdays and distances
        calculateGenerationValues(allAdults, cConstants->missionObjectives, objectiveAvgValues, duplicateNum, minDistance, avgDistance, maxDistance, minSteps, avgSteps, maxSteps, generation, avgAge, avgBirthday, oldestBirthday);

        //Assumes oldAdults is in rankDistance order
        changeAnneal (oldAdults, cConstants, currentAnneal, generation);

        //Get the run's total time
        totRunTime = (std::chrono::system_clock::now() - runStartTime);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE RECORD-_-_-_-_-_-_-_-_-_\n\n";
        //Print out necessary info for this generation
        genOutputs.printGeneration(cConstants, allAdults, objectiveAvgValues, generation, currentAnneal, numErrors, duplicateNum, minSteps, avgSteps, maxSteps, minDistance, avgDistance, maxDistance, avgAge, generation-oldestBirthday, avgBirthday, oldestBirthday, (totRunTime.count()/(generation+1)));

        // Before replacing new adults, determine whether all are within tolerance
        // Determines when loop is finished
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CONVERGENCE CHECK-_-_-_-_-_-_-_-_-_\n\n";
        convergence = checkTolerance(oldAdults, cConstants);

        //std::sort(allAdults.begin(), allAdults.end(), LowerMarsDist);
        //if (allAdults[0].minMarsDist < MSOI*cConstants->MSOI_scale && worstOPD < allAdults[0].orbitPosDiff) {
        //    std::string tempPath = genOutputs.outputPath;
        //    genOutputs.outputPath = tempPath+"badAdult\\";
        //    mkdir(genOutputs.outputPath.c_str());
        //    genOutputs.finalRecord(cConstants, allAdults[0], generation);
        //    genOutputs.outputPath = tempPath;
        //    worstOPD = allAdults[0].orbitPosDiff;
        //    std::cout << "\nBAD ADULT PRINTED\n\n";
        //}
        std::sort(allAdults.begin(), allAdults.end(), rankRaritySort);        
        
        //Increment the generation counter
        ++generation;
    
        //Loop exits based on result of checkTolerance and if max_generations has been hit
    } while ( !convergence && generation < cConstants->max_generations);

    //Handle the final printing
    genOutputs.printFinalGen(cConstants, allAdults, gpuValues, convergence, generation, numErrors, duplicateNum, oldestBirthday, (totRunTime.count()/(generation+1))); 

    return calcPerS;
}