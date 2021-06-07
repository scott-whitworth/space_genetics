#include <iostream>  // cout
#include <iomanip> // setprecision(int)  

// Constructor used to initialize the earth calculation data
// Input: cConstants - Used to access impact date element data and the time range needed to be calculated (triptime max and min), passed into earthInitial_incremental()
EarthInfo::EarthInfo(const cudaConstants* cConstants) {
    // Setting up initial information
    startTime = cConstants->triptime_min; // Starting time (s), chronologically this is closest to impact time (0 would be exactly impact date)
    endTime = cConstants->triptime_max;   // Ending time (s), chronologically this is earliest time away from impact date
    timeRes = cConstants->timeRes;        // Time resolution for storing data points (s)
    tolData = ((endTime-startTime)/timeRes) + 1; // Total Number of Data points in earthCon based on duration in seconds divided by resolution, plus one for the last 'section'

    // Alocate memory for earthCon, one entry for every data point
    earthCon = new elements<double> [tolData];

    // Assigning the position of the earth at impact to variable earth. Passed into earthInitial_incremental and rk4Reverse.
    elements<double> earth = elements<double>(cConstants->r_fin_earth, cConstants->theta_fin_earth, cConstants->z_fin_earth, cConstants->vr_fin_earth, cConstants->vtheta_fin_earth, cConstants->vz_fin_earth);

    // Get the initial position and velocity of the earth from startTime away from impact date.
    earth = earthInitial_incremental(0, startTime, earth, cConstants);

    // Setting the first element of earthCon to be equal to the earth conditions at startTime away from impact
    earthCon[0] = earth;

    // Shows progress of earth position calculations before the optimization in cuda can occur.
    std::cout << "Calculating earth positions for the trip time range" << std::endl;
    std::cout << "          10 20 30 40 50 60 70 80 90 100" << std::endl;
    std::cout << "progress:[";

    // Iterate backwards until total number of data points acquired
    for (int i = 1; i < tolData; i++) { 
        // Calculates earth's condition at each point (time) from the previously calculated point.
        earthCon[i] = earthInitial_incremental(calc_time(i)-timeRes, calc_time(i), earth, cConstants); // Obtaining conditions of the earth

        // Filling progress bar, 30 value derived from the fact that there is 30 character spaces to fill in progress bar
        if ((i % (tolData/30)) == 0) {
            std::cout << ">";
        }

        // Sets earth equal to the conditions calculated for a given time to be used as reference point in next usage of earthInitial_incremental
        earth = earthCon[i];
    }
    // Closing progress bar and adding a couple of empty line spaces
    std::cout << "]\n\n";
}

// Returns conditions of earth for a given time input
// Input: currentTime - time offset from impact backwards in time (larger value refers further back) in units of seconds
// Output: Returns an element that is to earth's position/velocity at currentTime away from impact using interpolate
//         as currentTime very likely does not directly corelate to an explicit derived element in earthCon
elements<double> EarthInfo::getCondition(const double & currentTime) {
    // Defining variables to use
    elements<double> lower;
    elements<double> upper;
    elements<double> result;
    
    int index;
    // Setting index equal to the nearest index of data based on time
    index = calcIndex(currentTime); 
    
    if (index < calcIndex(startTime)) {
        std::cout << "Earth condition index being accessed is less than startTime! Returning nearest index\n";
        index = calcIndex(startTime);
    }
    else if (index > calcIndex(endTime) - 1) {
        std::cout << "Earth condition index being accessed is greater than to endTime-1! Returning nearest valid index\n";
        index = calcIndex(endTime) - 1; // is (index-1) because the index variable is the lower of two, the keeps from going out of bounds with the upper index value 
    }

    // The lower index weight is equal to the index
    lower = earthCon[index];
    // The upper index weight is equal to the index plus 1
    upper = earthCon[index + 1];

    // Weights are used for interpolation
    double lowerWeight = 1 - ((currentTime-calc_time(index)) / timeRes);
    double upperWeight = 1 - ((calc_time(index+1)-currentTime) / timeRes);

    // Derive the element of Earth at this time using interpolate
    result = interpolate(lower,upper,lowerWeight,upperWeight);

    return result;
}

// Takes in a time and outputs a corresponding index (location of data).
int EarthInfo::calcIndex(const double & currentTime) {
   return static_cast<int>((currentTime-startTime)/timeRes);
}

// Takes in an index and outputs the time corresponded to that index.
double EarthInfo::calc_time(const int & currentIndex) {
    return startTime + (currentIndex*timeRes);
}

// Returns the tolData
int EarthInfo::getTolData() {
    return tolData;
}

// Takes in the lower index, upper index, and their weights in order to calculate an approximate Earth element for a time between two index.
// Input: lower & upper - elements that the interpolate is between of
//        lowerWeight & upperWeight - weight values to adjust the interpolation in cases where the time is closer to one index than another (not midpoint)
elements<double> EarthInfo::interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight) {
    return (lower*lowerWeight)+(upper*upperWeight);
}

// Destructor, deallocates earthCon data
EarthInfo::~EarthInfo() {
    delete [] earthCon;
}

// Determines the next step back from a given element using rk4Reverse, used in the constructor of EarthInfo
// Input for earthInitial_incremental:
//      timeInitial: time (in seconds) of the impact date
//      tripTime: optimized time period of the overall trip
//      earth: earth's position and velocity on the impact date October 5, 2022 (passed in from earthInfo.cpp)
elements<double> earthInitial_incremental(double timeInitial, double tripTime, const elements<double> & earth, const cudaConstants * cConstants) {
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = (tripTime - timeInitial)/static_cast <double> (60); 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial, tripTime, earth, deltaT, yp, cConstants->rk_tol, cConstants);
 
  return yp;
}
