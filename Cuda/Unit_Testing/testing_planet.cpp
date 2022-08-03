#define UNITTEST
#include "../Planet_calculations/planetInfo.h"  // For launchCon and EarthInfo()
#include <math.h>  //allows us to access isnan
#include <random> //for rng
#include <time.h>

//NOTE: Because there were so few of these, we did not bother to make a .h for this,
//      but if more planet unit testing is done, it would like be worth it to move these
//      function headers to a .h file
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

    utcConstants->orbitalPeriod = 1.0e8; // so it is easy to see the period and it should be the same length as the planet's position is (s) 

    //ALL THE STUFF SO EARTH'S INITIAL POSITION CAN BE CALCULATED
    utcConstants->triptime_max=2.0*SECONDS_IN_YEAR; //set to 2 years
    utcConstants->triptime_min=1.0*SECONDS_IN_YEAR; //set to 1 year
    utcConstants->timeRes=1800; // Earth Calculations Time Resolution Value
    
    //Bennu config values but shortened to make comparisons easier
    utcConstants->r_fin_earth=9.86E-01; //in AU
    utcConstants->theta_fin_earth=1.25; //in rad
    utcConstants->z_fin_earth=-4.33E-05; //in AU
    utcConstants->vr_fin_earth=-1.8E-09; //in AU/s
    utcConstants->vtheta_fin_earth=2.0E-07; //in AU/s
    utcConstants->vz_fin_earth=-2.3E-12; // in AU/s

    //Bennu Mars values
    utcConstants->r_fin_mars=1.422486429200229; //in AU
    utcConstants->theta_fin_mars=0.433301411341479; //in rad
    utcConstants->z_fin_mars=-1.916394400307147E-02; //in AU
    utcConstants->vr_fin_mars=1.141483741401610E-08; //in AU/s
    utcConstants->vtheta_fin_mars=1.719277837221393E-07; //in AU/s
    utcConstants->vz_fin_mars=4.887178957301878E-09; //in AU/s
    

// SETTING A RANDOM NUMBER GENERATOR (rng) TO BE USED BY FUNCTIONS
    // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    std::mt19937_64 rng(utcConstants->time_seed);

    //VERY complicated part of the code with some possibility of errors -> just needed for the child constructor with rkParameters and cConstants as its arguments
    marsLaunchCon = new PlanetInfo(utcConstants, MARS); 

// CALLING THE DIFFERENT UNIT TESTING ALGORITHMS
    bool allWorking = true;

    if (testGetConditions(utcConstants)){
        std::cout << "Getting the position of Mars works as expected - the part inside the parantheses is how many seconds before the end of the mission" << std::endl;
    }
    else{
        std::cout << "Something about getting the position of Mars did not work as expected" << std::endl;
    }

    delete utcConstants;
    delete marsLaunchCon;
    return allWorking;
}

bool testGetConditions(const cudaConstants* utcConstants){
    bool expectedResults = true;

    //uses get conditions to try and get the final conditions of Mars (the values from the config)
    //since time here goes backwards (starting from the asteroid and counting back in time), 0 is at the asteroid 
    //while something like 100 would mean 100 seconds before the spacecraft reached the asteroid
    elements<double> final = (*marsLaunchCon).getCondition(0);
    
    //uses the overloaded << to print the position of Mars at this time to the terminal
    std::cout << final << std::endl;

    //this returns whether or not we're getting an expected result
    //Mars's position should be its final position
    if (!(final.r == utcConstants->r_fin_mars)){

        expectedResults = false;
    }

    //this should be the position of Mars an hour before the spacecraft reaches its target
    elements<double> hourBefore = (*marsLaunchCon).getCondition(3600);

    //Positions of Mars on 2018-Dec-03 16:00:00.0000 according to https://ssd.jpl.nasa.gov/horizons/app.html#/
    double rMarsHourBefore = 1.4224453435;
    double thetaMarsHourBefore = 0.4328662874;
    double zMarsHourBefore = -1.9181536137E-02;


    //uses the overloaded << to print the position of Mars at this time to the terminal
    std::cout << hourBefore << std::endl;

    //this returns whether or not we're getting an expected result
    //Mars's position should be its final position
    if (!(hourBefore.r < rMarsHourBefore+1E-10 && hourBefore.r > rMarsHourBefore-1E-10)){
        std::cout << "The radius is wrong when looking at an hour before the end of the mission... " << std::endl;
        expectedResults = false;
    }
    if (!(hourBefore.theta < thetaMarsHourBefore+1E-10 && hourBefore.theta > thetaMarsHourBefore-1E-10)){
        std::cout << "The theta is wrong when looking at an hour before the end of the mission... " << std::endl;
        expectedResults = false;
    }
    if (!(hourBefore.z < zMarsHourBefore+1E-12 && hourBefore.z > zMarsHourBefore-1E-12)){
        std::cout << "The z value is wrong when looking at an hour before the end of the mission... " << std::endl;
        expectedResults = false;
    }

    double tripTime = 47304000; //Makes a fake 1.5 year trip time (don't use 1.5*SECONS_IN_YEAR because that is abou 9 hours off)

    //this should be the position of Mars an when the spacecraft launches (1.5 years before the spacecraft reaches its target)
    elements<double> launchTime = (*marsLaunchCon).getCondition(tripTime);

    //uses the overloaded << to print the position of Mars at this time to the terminal
    std::cout << launchTime << std::endl;

    //Mars coordinates from 2017-Jun-03 17:00:00.0000 (1.5 years before) 
    double rMarsLaunch = 1.5912271219;
    double thetaMarsLaunch = 1.7327774331;
    double zMarsLaunch = 3.9206489815E-02;

    if (!(launchTime.r < rMarsLaunch+1E-10 && launchTime.r > rMarsLaunch-1E-10)){
        std::cout << "The radius is wrong when looking at the spacecraft's launch... " << std::endl;
        std::cout << "The radius is off by " << launchTime.r - rMarsLaunch << " AU" << std::endl;
        expectedResults = false;
    }
    if (!(launchTime.theta < thetaMarsLaunch+1E-10 && launchTime.theta > thetaMarsLaunch-1E-10)){
        std::cout << "The theta is wrong when looking at the spacecraft's launch... " << std::endl;
        std::cout << "Theta is off by " << launchTime.theta - thetaMarsLaunch << " rad" << std::endl;
        expectedResults = false;
    }
    if (!(launchTime.z < zMarsLaunch+1E-12 && launchTime.z > zMarsLaunch-1E-12)){
        std::cout << "The z value is wrong when looking at the spacecraft's launch... " << std::endl;
        std::cout << "z is off by " << launchTime.z - zMarsLaunch << " AU" << std::endl;
        expectedResults = false;
    }

    elements<double> usingConditionDev = getConditionDev(tripTime, utcConstants, marsLaunchCon);

    //uses the overloaded << to print the position of Mars at this time to the terminal
    std::cout << usingConditionDev << std::endl;

    if (launchTime.r != usingConditionDev.r){
        std::cout << "using getConditionDev is not getting the same answer as getCondition" << std::endl;
        std::cout << "The difference in their radii is " << launchTime.r - usingConditionDev.r << std::endl;
        expectedResults = false;
    }

    return expectedResults;
}