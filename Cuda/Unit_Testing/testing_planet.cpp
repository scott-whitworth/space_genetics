#define UNITTEST
#include "../Planet_calculations/planetInfo.h"  // For launchCon and EarthInfo()
#include <math.h>  //allows us to access isnan
#include <random> //for rng
#include <time.h>

//runs the unit tests for the planets functions
bool runPlanetTests(bool printThings);

//Prints the conditions of the Earth at triptime_min, which should be the start point
bool testGetConditions(const cudaConstants* utcConstants);

bool runPlanetTests(bool printThings){
// SETTING UP CUDA CONSTANTS TO BE USED BY OTHER FUNCTIONS 
    //making cudaConstants to control what is going on in genetic algorithm while still using the original functions
    cudaConstants* utcConstants = new cudaConstants(); 

    // Seed used for randomization rng things, using seed 0 for consistancy / tracability 
    utcConstants->time_seed = 0; 

    //Flag for what type of landing you want. Current modes: "soft"=2, "hard"=1
    //testing for soft because as of 2022 that is the main mission type we are exploring 
    //additionally, mission type should have little impact on the genetics algorithms
    utcConstants->missionType = 2; 

    utcConstants->orbitalPeriod = 1.0e8; // so it is easy to see the period and it should be the same length as the planet's position is (s) 

    //ALL THE STUFF SO EARTH'S INITIAL POSITION CAN BE CALCULATED
    utcConstants->triptime_max=2.0e8; //set to 200 million seconds to make seeing this output easier on myself
    utcConstants->triptime_min=1.0e8; //set to 100 million seconds to make seeing this output easier on myself
    utcConstants->timeRes=1000; // Earth Calculations Time Resolution Value
    
    //Bennu config values but shortened to make comparisons easier
    utcConstants->r_fin_earth=9.86E-01;
    utcConstants->theta_fin_earth=1.25;
    utcConstants->z_fin_earth=-4.33E-05;
    utcConstants->vr_fin_earth=-1.8E-09;
    utcConstants->vtheta_fin_earth=2.0E-07;
    utcConstants->vz_fin_earth=-2.3E-12;

// SETTING A RANDOM NUMBER GENERATOR (rng) TO BE USED BY FUNCTIONS
    // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    std::mt19937_64 rng(utcConstants->time_seed);

    //VERY complicated part of the code with some possibility of errors -> just needed for the child constructor with rkParameters and cConstants as its arguments
    launchCon = new PlanetInfo(utcConstants, EARTH); 

// CALLING THE DIFFERENT UNIT TESTING ALGORITHMS
    bool allWorking = true;

    if (testGetConditions(utcConstants)){
        std::cout << "Good it works as expected" << std::endl;
    }
    else{
        std::cout << "IT DID NOT WORK" << std::endl;
    }

    delete utcConstants;
    delete launchCon;
    return allWorking;
}

bool testGetConditions(const cudaConstants* utcConstants){
    //uses get conditions to try and get the initial conditions of the Earth
    elements<double> initial = (*launchCon).getCondition(utcConstants->triptime_max-1);
    
    //uses the overloaded << to print this to the terminal
    std::cout << initial << std::endl;

    //this returns whether or not we're getting an expected result
    if (initial.r == utcConstants->r_fin_earth){
        return true;
    }

    return false;
}