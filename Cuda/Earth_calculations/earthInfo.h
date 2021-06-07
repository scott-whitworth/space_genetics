#ifndef earthInfo_h
#define earthInfo_h

#include "../Motion_Eqns/elements.h"
#include "../Config_Constants/config.h"
#include "../Runge_Kutta/runge_kutta.h"

class EarthInfo {
    private:
        // Elements pointer array that holds all of the earth conditions for a given time range
        elements<double> *earthCon;
        // First time point for a given time span offseted from impact date, units of seconds
        double startTime;
        // Last time point for a given time span offseted from impact date, units of seconds
        double endTime;
        // The resolution of data points (ex: 3600 = hours, 60 = minutes, 1 = seconds)
        double timeRes;
        // The total amount of data points for a run. Calculated by time span divided by timeRes. Is the length of earthCon
        int tolData;

        // Takes in a time and outputs a corresponding index (location of data).
        int calcIndex(const double & currentTime);
        // Takes in an index and outputs the time corresponded to that index.
        double calc_time(const int & currentIndex);

        // Takes in the lower index, upper index, and their weights in order to calculate an approximate Earth element for a time between two index.
        // Input: lower & upper - elements that the interpolate is between of
        //        lowerWeight & upperWeight - weight values to adjust the interpolation in cases where the time is closer to one index than another (not midpoint)
        elements<double> interpolate(const elements<double> & lower, const elements<double> & upper, const double & lowerWeight, const double & upperWeight);

    public:
        // Constructor used to initialize the earth calculation data
        // Input: cConstants - Used to access impact date element data and the time range needed to be calculated (triptime max and min), passed into earthInitial_incremental()
        EarthInfo(const cudaConstants* cConstants);
        
        // Returns the interpolated conditions of earth for a given time input, using interpolate if currentTime does not directly corelate to an explicit derived element in earthCon
        // Output: Returns an element that is to earth's position/velocity at currentTime away from impact (backwards)
        elements<double> getCondition(const double & currentTime);

        // Returns the total amount of data for a run with a given time span and resolution.
        int getTolData();

        // Clears dynamic memory, used at end of program when optimize run is completed
        // Deallocates earthCon data
        ~EarthInfo();
};

// Determines the next step back from a given element using rk4Reverse, used in the constructor of EarthInfo
// Input for earthInitial_incremental:
//      timeInitial: time (in seconds) of the impact date
//      tripTime: optimized time period of the overall trip
//      earth: earth's position and velocity on the impact date October 5, 2022 (passed in from earthInfo.cpp)
elements<double> earthInitial_incremental(double timeInitial, double tripTime,const elements<double> & earth, const cudaConstants * cConstants);

#include "earthInfo.cpp"

// Global variable for launchCon (assigned content in optimization.cu)
EarthInfo *launchCon;

#endif