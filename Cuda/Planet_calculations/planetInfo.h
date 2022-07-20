#ifndef planetInfo_h
#define planetInfo_h

#include "../Motion_Eqns/elements.h"
#include "../Config_Constants/config.h"
#include "../Runge_Kutta/runge_kutta.h"

class PlanetInfo{
    //TODO: Switch public designation back to private once we are finished running tests
    //private:
    public:
        // Elements pointer array that holds all of the planet conditions for a given time range
        // This array starts with the position of the planet when the spacecraft is at the asteroid 
        //      planetCon[0] = position of planet when the spacecraft is at the asteroid
        //      planetCon[tripTime] = position of the planet when the spacecraft is about to leave Earth and begin its trip
        //      planetCon[T] = position of the planet T timesteps (of cConstants->timeRes seconds) before the spacecraft reaches the asteroid
        // Thus every time step technically goes backwards in time -> the next element of the array is 1 hour BEFORE 
        //      the one that follows it. 
        // This is the case because the actual calendar date when the spacecraft lands on the asteroid is AFTER
        //      the calendar date of when the spacecraft leaves the Earth, but in this array the position of the 
        //      planet when the spacecraft is at the asteroid comes FIRST
        // This is can be kind of confusing on a mission starting from Earth and heading to the asteroid
        elements<double> *planetCon;
        //enum planet status that is given by the constuctor when it is called specifying which planet is being calculated
        int planetName;
        // First time point for a given time span offseted from impact date, units of seconds
        double startTime;
        // Last time point for a given time span offseted from impact date, units of seconds
        double endTime;
        // The resolution of data points (ex: 3600 = hours, 60 = minutes, 1 = seconds)
        double timeRes;
        // The total amount of data points for a run. Calculated by time span divided by timeRes. Is the length of planetCon
        int tolData;

        // Takes in a time and outputs a corresponding index (location of data).
        int calcIndex(const double & currentTime);
        // Takes in an index and outputs the time corresponded to that index.
        double calc_time(const int & currentIndex);

        // Takes in the lower index, upper index, and their weights in order to calculate an approximate Planet element for a time between two index.
        // Input: lower & upper - elements that the interpolate is between of
        //        lowerWeight & upperWeight - weight values to adjust the interpolation in cases where the time is closer to one index than another (not midpoint)
        elements<double> interpolate(const elements<double> & lower, const elements<double> & upper, const double & lowerWeight, const double & upperWeight);

        // Uses a planetStatus to determine which planet's position at the time when the spacecraft reaches the asteroid we want planet to be filled with
        // Input: cConstants - to access the final positions and velocities of the planet as specified in the config file
        //        planet - a set of elements that will be filled with the planet's position when the spacecraft reaches the asteroid
        //        planetStatus - allows us to identify which planet we are  
        void getPlanetPosition(const cudaConstants* cConstants, elements<double> & planet, const int & planetStatus);

    public:

        // Constructor used to initialize the planet calculation data
        // Input: cConstants - Used to access impact date element data and the time range needed to be calculated (triptime max and min), passed into planetInitial_incremental()
        //        planetStatus - tells us which planet we are initializing because with Earth, for example, we only need to know its positions in the launch window
        //                      not for the rest of the trip
        PlanetInfo(const cudaConstants* cConstants, const int & planetStatus);
        
        // Returns the interpolated conditions of a planet for a given time input, using interpolate if currentTime does not directly corelate to an explicit derived element in planetCon
        // Input: currentTime - seconds BEFORE the spacecraft reaches the asteroid
        // Output: Returns an element that is to a planet's position/velocity at currentTime away from impact (backwards)
        // IMPORTANT NOTE: 
        //      currentTime is seconds BEFORE the spacecraft reaches the ASTEROID
        //      currentTime is NOT the number of seconds that have passed since the spacecraft left Earth 
        // NOTE: There is a device version of getCondition -> if changes are made here, you will likely want to make them there too 
        elements<double> getCondition(const double & currentTime);
        

        // Returns the total amount of data for a run with a given time span and resolution.
        // In other words, it returns the size of planetCon
        int getTolData();

        // Returns a pointer to the planetCon array 
        // Allows you to move the positions of your planets at different time steps through the code without incurring 
        //      the memory cost of all the other elements of planetInfo
        elements<double>* getAllPositions();

        // Clears dynamic memory, used at end of program when optimize run is completed
        // Deallocates planetCon data
        ~PlanetInfo();

};


// Reports size/number of elements in the planetCon array -> passed into Malloc
// Input: cConstants - to access triptime_max and timeRes
// Output: 6 times the number of elements in the planetCon array 
//         This number is multiplied by 6 because that is how many things are in elements<T>
int getPlanetSize(const cudaConstants * cConstants);

// Determines the next step back from a given element using rk4Reverse, used in the constructor of PlanetInfo
// Input for planetInitial_incremental:
//      timeInitial: time (in seconds) of the impact date
//      tripTime: optimized time period of the overall trip
//      planet: planet's position and velocity on the arrival date (passed in from planetInfo.cpp)
elements<double> planetInitial_incremental(double timeInitial, double tripTime,const elements<double> & planet, const cudaConstants * cConstants);

// Interpolates the data to find the position of a planet currentTime seconds BEFORE the spacecraft reaches the asteroid
// IMPORTANT NOTE: 
//      currentTime is seconds BEFORE the spacecraft reaches the ASTEROID
//      currentTime is NOT the number of seconds that have passed since the spacecraft left Earth
// Input: currentTime - seconds before the spacecraft reaches the asteroid
//        cConstants - to access triptime_max and timeRes
//        planetConditions - the array of elements that holds the position of the planet at timeRes intervals
// Output: The position of the planet currentTime seconds BEFORE the spacecraft reaches the asteroid
// NOTE: This is a device version of getCondition -> if changes are made to getCondition, you will likely want to make them here too 
__host__ __device__ elements<double> getConditionDev(const double & currentTime, const cudaConstants * cConstants, const elements<double>* planetConditions);


#include "planetInfo.cpp"
// Global variable for launchCon (assigned content in optimization.cu)
PlanetInfo *launchCon;
PlanetInfo* marsLaunchCon;

#endif