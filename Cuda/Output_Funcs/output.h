#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream> //for file input/output
#include <vector> //allows for using vectors instead of just dynamic arrays

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// input: Individual to be displayed (assumed to be the best individual of the pool) 
//        and the value for the current generation iterated
// output: display: generation, individual's posDiff, speedDiff
void terminalDisplay(const Adult& individual, unsigned int currentGeneration);

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
void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, double *work, double *dE, double *Etot_avg, const cudaConstants* config);

// Main Output, final results of genetic algorithm
// input: x[] - array of OPTIM_VARS for a single individual
//        generation - generation num of individual        
//        cConstants - Access constants info such as target element, earth element, derive spaceCraft element, also other values such as rk_tol
//        best - To access the best individual (pool[0])
// output: file orbitalMotion-[time_seed].bin is created that holds spacecraft RK steps and error
//         file finalOptimization-[time_seed].bin is created that holds earth/ast/ and trajectory parameter values
void trajectoryPrint(double x[], int generation, const cudaConstants* cConstants, const Adult best);

// Record progress of individual
// Called if record_mode is true at end of optimize process
// input: generation - the positional performance of the individual
//        numStep - the number of steps taken in the RK method
//        start - the array of optimized initial parameters, OPTIM_VARS size
//        yp - the final RK solution at numStep
//        config - cudaConstants object for accessing thruster_type information
// output: output file (progressiveAnalysis.csv) is appended individual values/parameter information
void progressiveAnalysis(int generation, int numStep, double *start, elements<double> & yp, const cudaConstants *config);

// Initialize some of the files used in record mode with header rows
// input: cConstants - to access time_seed for deriving file name conventions
// output: files genPerformanceT-[time_seed].csv is given initial header row info for generation, best posDiff, best speedDiff, and parameters of best individual
void initializeRecord(const cudaConstants * cConstants);

// Take in the current state of the generation and appends to files
// assumes initializeRecord() had already been called before (therefore no need to output a header row)
// input: cConstants - access time_seed to derive file name
//        pool - access individuals of current generation, assumed to be ordered by cost
//        generation - record current generation
//        new_anneal - outputs the current generation's anneal
//        anneal_min - outputs the current generation's anneal_min
//        poolSize - to know size of the pool, currently unused
// output: genPerformanceT-[time_seed].csv is appended parameter information on the best individual in pool
void recordGenerationPerformance(const cudaConstants * cConstants, const std::vector<Adult>& pool, double generation, double new_anneal, int poolSize, double anneal_min);

// Method for doing recording information at the end of the optimization process
// input: cConstants - to access config info
//        pool - To access the best individual (pool[0])
//        generation - to record the generation value 
// output: trajectoryPrint is called to generate binary files for MATLAB plotting
void finalRecord(const cudaConstants* cConstants, const std::vector<Adult> pool, int generation);

// ** Currently unused ** 
// input: cConstants - access time_seed to derive file name
// output: mutateFile[time_seed].csv is given a header row, now ready to be used for progressiveRecord()
void setMutateFile(const cudaConstants* cConstants);

// ** Currently not used ** 
// input: cConstants - access time_seed to derive file name
//        generation - for storing into the file, the generation that the mutation is occuring
//        annealing - for storing into the file, the anneal value being used
//        numGenes - the number of genes that are to be mutated
//        recordLog[OPTIM_VARS] - a "mutate mask" that is an array corresponding to the optimization array that contains the random mutate values being applied
// output: file mutateFile-[time_seed].csv is appended a new row containing mutate information
void recordMutateFile(const cudaConstants * cConstants, double generation, double annealing, int numGenes, double recordLog[OPTIM_VARS]);

// method that stores information of launchCon of timeRes*24 resolution
// input: cConstants - access time range and resolution info
//        launchCon (global variable) - access elements of earth 
// output: EarthCheckValues.csv is created and holds rows of element info on earth with timeStamp on each row
void recordEarthData(const cudaConstants * cConstants);

// Takes in a pool and records the parameter info on all individuals
// input: cConstants - to access time_seed in deriving file name
//        pool - holds all the individuals to be stored
//        poolSize - to use in iterating through the pool
//        generation - used in deriving file name
// output: file generation#[generation]-[time_seed].csv is created with each row holding parameter values of individuals
void recordAllIndividuals(std::string name, const cudaConstants * cConstants, const std::vector<Adult>& pool, int poolSize, int generation);

// Record initial and final fuel masses along with tripTime and relative velocity at impact
// input: cConstants - access config constants
//        solution - best individual parameters from the final pool
//        fuelSpent - total fuel spent
//        best - To access the best individual (pool[0])
// output: fuelOutput.csv - output file holding fuel consumption and impact data
void recordFuelOutput(const cudaConstants* cConstants, double solution[], double fuelSpent, const Adult best);


#include "output.cpp"

#endif
