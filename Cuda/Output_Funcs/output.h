#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream> //for file input/output
#include <vector> //allows for using vectors instead of just dynamic arrays

//Structure will hold the run's output folder address, many output functions, and will handle generational & final outputs
//This will help make output folders more flexible
struct output 
{
    //String stores the overall folder where each individual run's files are stored
    //Will be ..\Output_Files\ unless set otherwise
    std::string baseFolder;

    //Stores the folder address for the individual run's output files
    //Will be derived from baseFolder
    std::string outputPath; 

    //Base constructor, should never be used
    output():baseFolder("..\\Optimization\\") {} 

    //Constructor will serve as both a default constructor and a specific constructor
    //  It will allow for a specific output base to be specified, but will assume a standard base folder if none are specified
    //Inputs: cConstants - program cudaConstants
    //        selectFolder - the desired folder to hold the results of the runs, defaults to Output_Files
    //Output: The folders & files necessary for this run will be prepared
    output(const cudaConstants* cConstants, const std::string& selectFolder = "..\\Output_Files\\"); 

    // Function will handle the generation-to-generation outputs
    // Inputs: Passed to printBestAdults(), recordGenerationPerformance(), recordGenSimple(), and recordAllIndividuals()
    //              See those functions' header files for detail on how the inputs are used
    // Outputs: The functions lised above will be called depending on the current generation and cConstants's write frequency variables
    void printGeneration(const cudaConstants * cConstants, const std::vector<Adult>& allAdults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, int& errorNum, const int& duplicateNum, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday); 

    // Function will handle printing at the end of a run
    // Inputs: Passed to recordAllIndividuals(), printBestAdults(), reportRun(), & finalRecord()
    //              See those functions' header files for detail on how the inputs are used
    // Outputs: printBestAdults() and reportRun() will be called regardless, with recordAllIndividuals() and finalRecord() being called conditionally
    void printFinalGen(const cudaConstants * cConstants, std::vector<Adult>& allAdults, const bool& converged, const int& generation, int& errorNum, const int& duplicateNum, const int& oldestBirthday, const float& avgGenTime); 
    
    // Initialize genPerformance with header rows
    // input: cConstants - to access time_seed for deriving file name conventions
    // output: files genPerformanceT-[time_seed].csv is given initial header row info for generation, staring parameters, outputs, and age values for the best adults along with generation-wide values
    void initializeGenPerformance(const cudaConstants * cConstants);

    // Initialize simpleGenPerformance with header rows
    //      This file will operate as a similified version of genPerformance for quick progress tracking
    // Inputs: cConstants - the cudaConstants, used to access objectives
    // Output: the simpleGenPerformance + seed will be created, with only generation, parameter values, and progress for the best individual
    void initializeSimpleGenPerformance(const cudaConstants* cConstants);
    
    // Checks if needed output files exist, creates them if they don't exist, and copies the run's config files to them
    // Inputs: cConstants - the cuda constants, will be used to access the destination
    // Output: the baseFolder folder with a sub-folder named after the time seed will be created, with the sub-folder being filled with the config files
    void prepareFolders(const cudaConstants* cConstants);

    // Assist function to prepareFolders, which will copy the config files into a seed output folder
    // Inputs: source - the source config file
    //         desination - the destination of the source file
    //              NOTE: must include the desired file name for the destination as well, it can't just end at a folder, or the source won't be copied
    // Output: The source file will be put in the desination
    void copyFile (const std::string& source, const std::string& destination);

    // Take in the current state of the generation and appends to files
    // assumes initializeGenPerformance() had already been called before (therefore no need to output a header row)
    // inputs: cConstants - access time_seed to derive file name
    //         adults - A vector of adults which is used to pull the best adults from; not passed by reference to maintain the rank-distance sort of the original vector
    //         objectiveAvgValues - a vectors with a size matching the number of objectives
    //               Each index holds the average generation value of the parameter associated with the matching objective
    //         generation - record current generation
    //         new_anneal - outputs the current generation's anneal
    //         errorNum - the amount of errors found this generation
    //         duplicateNum - the number of duplicate adults found
    //         minDist - minimum distance for the generation
    //         avgDist - average distance for the generation
    //         maxDist - the max distance for the generation
    //         avgAge - the generation-adjusted average generation age
    //         oldestAge - the generation-adjusted oldest generation age
    //         avgBirthday - the non-generation-adjusted average generation age
    //         oldestBirthday - the non-generation-adjusted oldest generation age
    // output: genPerformanceT-[time_seed].csv is appended parameter information on the best individual in pool
    void recordGenerationPerformance(const cudaConstants * cConstants, std::vector<Adult> adults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, const int& errorNum, const int& duplicateNum, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday);

    // Record highlights of the full genPerformance file
    //  Note: assumes initializeRecord() has been called
    // Inputs: cConstants - the cudaConstants, used to access objectives
    //         adults - A vector of adults which is used to pull the best adults from; not passed by reference to maintain the rank-distance sort of the original vector
    // Output: The essential stats for the generation's best adult is added to simpleGenPerformance
    void recordGenSimple (const cudaConstants* cConstants, std::vector<Adult> adults, const std::vector<double>& objectiveAvgValues, const int& generation);

    // Takes in a vector of adults and records the parameter info on all individuals
    // input: name - the desired name of the output, will be inserted into the file name
    //        cConstants - to access time_seed in deriving file name
    //        adults - vector of adults that will be recorded
    //        generation - used in deriving file name
    // output: file generation#[generation]-[time_seed].csv is created with each row holding parameter values of individuals
    void recordAllIndividuals(std::string name, const cudaConstants * cConstants, const std::vector<Adult>& adults, const int& generation);

    // Method for doing recording information at the end of the optimization process
    // input: cConstants - to access config info
    //        bestAdult - The final best individual
    //        generation - to record the generation value 
    // output: trajectoryPrint is called to generate binary files for MATLAB plotting
    void finalRecord(const cudaConstants* cConstants, const Adult& bestAdult, const int& generation);

    // Function which prints the general info of a run when its over
    // Inputs: cConstants - cuda constants, used for accessing seed and objectives
    //         oldAdults - the oldAdult vector, will report the parameters from the top adults within the vector
    //              NOTE: assumes oldAdults is in rank-distance sort
    //         converged - a bool which indicates whether the run converged or not
    //         generation - the final generation
    //         avgGenTime - the average time (in seconds) it took a generation to complete during the run
    // Output: a file with at-a-glance info of runs is appended to within Output_Files
    void reportRun(const cudaConstants* cConstants, const std::vector<Adult>& oldAdults, const bool& converged, const int& generation, const float& avgGenTime); 

    // Main Output, final results of genetic algorithm
    // input: x[] - array of OPTIM_VARS for a single individual
    //        generation - generation num of individual        
    //        cConstants - Access constants info such as target element, earth element, derive spaceCraft element, also other values such as rk_tol
    //        best - To access the best individual (pool[0])
    // output: file orbitalMotion-[time_seed].bin is created that holds spacecraft RK steps and error
    //         file finalOptimization-[time_seed].bin is created that holds earth/ast/ and trajectory parameter values
    void trajectoryPrint(double x[], int generation, const cudaConstants* cConstants, const Adult& best);

    // Records error in energy conservation due to thrust calculations
    // For a single set of rkParameters
    // input: time - array, time steps
    //        yp - array, position/velocity
    //        gamma - array, in-plane thrust angles
    //        tau - array, out-of-plane thrust angles
    //        lastStep - number of steps taken in RK
    //        accel - array, magnitudes of acceleration due to thrust
    //        fuelSpent - array, aggregate amounts of fuel spent
    //        wetMass - initial mass of spacecraft and fuel
    //        config - cudaConstants object for accessing genetic.config information
    // output: work - array of work done by thruster between time steps (J)
    //         dE - array of changes in mechanical energy between time steps (J)
    //         Etot_avg - array of average mechanical energy between time steps (J)
    // Used to fill orbitalMotion{}.bin
    void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, double *work, double *dE, double *Etot_avg, const cudaConstants* config, elements<double> *mars);

    // method that stores information of launchCon of timeRes*24 resolution
    // input: cConstants - access time range and resolution info
    //        launchCon (global variable) - access elements of planet 
    // output: PlanetCheckValues.csv is created and holds rows of element info on planet with timeStamp on each row
    void recordMarsData(const cudaConstants * cConstants, const int & generation);
    
    // method that stores information of launchCon of timeRes*24 resolution
    // input: cConstants - access time range and resolution info
    //        launchCon (global variable) - access elements of planet 
    // output: PlanetCheckValues.csv is created and holds rows of element info on planet with timeStamp on each row
    void recordEarthData(const cudaConstants * cConstants, const int & generation);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Below are functions associated with outputs, but could be run individually/without the output class
///////////////////////////////////////////////////////////////////////////////////////////////////////

// Function which will print the best rank-distance adult and the best adult for each objective to the terminal
// Inputs: cConstants - the cudaConstants, used for getting objectives
//         adults - the vector of adults which the best individuals will be pulled from; this function assumes no pre-done sort
//         generation - the current generation
//         numErrors - the number of errors in this generation
//         numDuplicates - the number of duplicates found in this generation
//         oldestBirthday - the oldest borthday from this generation
// Output: the top individual for rank distance and for each objective will have their stats printed to the terminal
void printBestAdults(const cudaConstants* cConstants, std::vector<Adult> adults, const int& generation, int& numErrors, const int& numDuplicates, const int& oldestBirthday); 

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// input:  Individual to be displayed (assumed to be the best individual of the pool) 
//         objectives - the program's current list of objectives
// output: display each objective's parameters for the passed in individual
void terminalDisplay(const Adult& individual, const std::vector<objective> objectives);

#include "output.cpp"

#endif
