#ifndef CONFIG_H
#define CONFIG_H

#include <iostream> // For << operator, fileRead, and std::cout
#include <random>
#include <string>
#include <vector>

#include "objective.h" // for storing objectives

// Structure that holds constant values related/used for the genetic algorithm that can be configured within a file 
// (as of August 5th 2020, that file address is hardcoded in main to be "../Config_Constants/genetic.config")
// Specifically allows for changing paramters without re-compiling
struct cudaConstants {
    double time_seed;    // Seed used for randomization within optimize function, if it's set to NONE the seed is set to time(0)
    int max_generations; // Maximum number of generations to evaluate in the genetic algorithm if not reaching a solution
    int run_count;       // How many runs of the optimization algorithm to perform using incremental seeds to not just repeat the same value (the actual change in the seed occurs in main)
    bool random_start;   // If set to false, the initial generation has individuals initialized from a file instead of randomly generated within a range
    bool record_mode;    // If set to true, functions that record information onto files such as genPerformance.csv.  The code still records a valid solution regardless of this setting
    std::string initial_start_file_address; // If random_start is false, use file_address to find what file is being used for the initial start

    //These variables are used in orbit missions to determine the desired final orbital radius and speed
    //They are initially set to -1 to allow the program to quickly determine if the goal is an orbit vs an impact/rendezvous
    double orbitalRadius = -1; // the radius of the orbit around a body
    double orbitalSpeed = -1; // the 

    int write_freq;       // Determine how many generations between calling recordGenerationPerformance() method (defined in Output_Funcs/output.cpp)
    int all_write_freq;   // Determine how many generations between calling recordAllIndividuals() method (defined in Output_Funcs/output.cpp)
    int disp_freq;        // Determine how many generations between calling terminalDisplay() method (defined in Output_Funcs/output.cpp)

    int best_count;        // Number of individuals that needs to be within the acceptable condition before ending the algorithm, also how many of the top individuals are recorded
    double anneal_initial; // initial value for annealing, meant to replace the previously used calculation involving ANNEAL_MIN and ANNEAL_MAX with something more simple
    double anneal_final;   // final value for annealing, anneal cannot be reduced beyond this point

    double mutation_amplitude; // The percentage for probability of mutating a gene in a new individual, called iteratively to mutate more genes until the check fails

    // Scalars used to modify the mutate_scales below, used to assist in making adults mutate more if needed
    // A value of 1 will have an individual's parameters mutate at the scale of the variables below
    double default_mutation_chance; 
    
    // Used in mutate(), affects the scale of change for the respective parameter values, in conjunction with annealing
    // Represents max bounds of mutation, mutation will never be +/- this value
    double gamma_mutate_scale; 
    double tau_mutate_scale; 
    double coast_mutate_scale;
    double triptime_mutate_scale;
    double zeta_mutate_scale;
    double beta_mutate_scale;
    double alpha_mutate_scale;

    // Used when random_start==true, is the max magnitude of the random parameter initial values (ranges from -value to +value)
    double gamma_random_start_range; // Sets all coefficients within the range
    double tau_random_start_range;
    double coast_random_start_range;
    double triptime_max; // Explicit bounds for valid triptime, impacts not just an individual's parameters but also range of EarthInfo's calculations regardless of random_start setting
    double triptime_min; // Explicit bounds for valid triptime, impacts not just an individual's parameters but also range of EarthInfo's calculations regardless of random_start setting
    double alpha_random_start_range;
    double beta_random_start_range; // For beta, only positive side of the range is used (0 to the value assigned)
    double zeta_random_start_range;

    // Used in thruster construction and corresponding calculations
    int thruster_type; // 0 is for no thruster, 1 is for NEXT ion thruster
    double dry_mass;   // Mass of the spacecraft with no fuel (kg)
    double fuel_mass;  // The mass quantity of fuel the spacecraft starts with, used to derive wet_mass
    double wet_mass;   // (derived) Wet mass is the total mass of the spacecraft (dry mass plus fuel), thinking make this derived from fuel_mass that would be in the config (kg)

    double coast_threshold; // 0 results in the thruster never coasting, 1 results in always coasting
    double c3scale;         // scalar multiplier for c3energy in constructor
    double c3energy;        // (initially assigned, but then multiplied with c3scale) specific energy of spacecraft at earth escape (m^2/s^2), determines vEscape
    double v_escape;        // (derived) magnitude of velocity at earth escape  (au/s), this variable is not directly accessible in the config file as it is derived from c3energy

    // The final position and velocity of the asteroid/target at impact date
    // Should be pulled from NASA database
    double r_fin_target;      // AU
    double theta_fin_target;  // Radians
    double z_fin_target;      // AU
    double vr_fin_target;     // AU/s
    double vtheta_fin_target; // AU/s
    double vz_fin_target;     // AU/s

    // The final position and velocity of the earth at impact date to be used as reference point
    double r_fin_earth;      // AU
    double theta_fin_earth;  // Radians
    double z_fin_earth;      // AU
    double vr_fin_earth;     // AU/s
    double vtheta_fin_earth; // AU/s
    double vz_fin_earth;     // AU/s

    // The final position and velocity of Mars at impact date to be used as reference point
    double r_fin_mars;      // AU
    double theta_fin_mars;  // Radians
    double z_fin_mars;      // AU
    double vr_fin_mars;     // AU/s
    double vtheta_fin_mars; // AU/s
    double vz_fin_mars;     // AU/s

    double v_impact; // AU/s, the official DART mission data, used in cost function of individuals to sort individuals with posDiff < pos_threshold

    double rk_tol;       // The relative/absolute (not sure which one it is) tolerance for the runge kutta algorithm
    double doublePrecThresh; // The smallest allowed double value for runge-kutta
    int cpu_numsteps;    // Constant higher precision for final CPU RK method, set to be equal to max_numsteps
    int min_numsteps;    // Minimum number of steps in runge kutta 
    int max_numsteps;    // Maximum number of steps in runge kutta 
    int num_individuals; // Number of individuals in the pool, each individual contains its own thread
    int survivor_count;  // Number of survivors selected, every pair of survivors creates 8 new individuals
    int thread_block_size;

    // Used in generating time range for Earth calculations (units in seconds), distance between explicit time intervals stored
    int timeRes;

    // minimum distance the spacecraft can be to the sun.
    double sun_r_min;
    
    //destination asteroid file. The constants for the asteroid/target and earth data are passed in with a different file in order to switch between asteroids/targets easier.
    std::string destination;

    //asteroid/target orbital constant
    double orbitalPeriod;

    double gravAssistDist;

    //Vector of objectives that will be solved by the program
    //  For more information of what the objective class holds, look at it's header file
    std::vector<objective> missionObjectives; 

    // Default constructor, assumes the files being pulled are genetic.config and mission.config
    cudaConstants();

    // Constructor, accepts a string argument for the config file path
    cudaConstants(std::string configFile);

    // Sets properties to what is within the config file
    // Input: File address that is used to open a text-based file and parses through to assign variables
    // Output: Properties explicitly set in the config file are set to values following equal sign, ignores comments or empty lines in files 
    // Notice: This does not verify much, if anything, about the config file!
    void FileRead(std::string fileName);

    // Transfers information from a config file line to the mission objectives vector
    //  Input:  a line taken from the mission config file that holds information on a desired mission objective
    //  Output: a new objective is pushed to the end of the missionObjectives vector
    void importObjective(std::string line);
};

// Output function to stream, with some formatting to help be more legible on terminal
std::ostream& operator<<(std::ostream& os, const cudaConstants& object);


#include "config.cpp"

#endif
