#include <iostream>  // cout
#include <iomanip> // setprecision(int)  

// Constructor used to initialize the planet calculation data
// Input: cConstants - Used to access impact date element data and the time range needed to be calculated (triptime max and min), passed into planetInitial_incremental()
PlanetInfo::PlanetInfo(const cudaConstants* cConstants, const int & planetStatus) {
    // Setting up initial information
    planetName = planetStatus;
    //earth starts at triptime as it will be a collection of impact postions
    //mars starts at zero as we care about all of its positions leading up to impact 
    if(planetStatus == EARTH){
        startTime = cConstants->triptime_min; // Starting time (s), chronologically this is closest to impact time (0 would be exactly impact date)
        endTime = cConstants->triptime_max;   // Ending time (s), chronologically this is earliest time away from impact date
    }else if(planetStatus == MARS){
        startTime = 0.0; // Starting time (s), chronologically this is closest to impact time (0 would be exactly impact date)
        endTime = cConstants->triptime_max;   // Ending time (s), chronologically this is earliest time away from impact date
    }
    
    timeRes = cConstants->timeRes;        // Time resolution for storing data points (s)
    tolData = ((endTime-startTime)/timeRes) + 1; // Total Number of Data points in planetCon based on duration in seconds divided by resolution, plus one for the last 'section'

    // Alocate memory for planetCon, one entry for every data point
    planetCon = new elements<double> [tolData];
    
    // Depending on the which planet this is, assigns the appropriate position of the planet at the spacecraft's arrival to planet
    elements<double> planet;
    getPlanetPosition(cConstants, planet, planetName);

    // Assigning the position of the earth at impact to variable earth. Passed into earthInitial_incremental and rk4Reverse.
    //elements<double> earth = elements<double>(cConstants->r_fin_earth, cConstants->theta_fin_earth, cConstants->z_fin_earth, cConstants->vr_fin_earth, cConstants->vtheta_fin_earth, cConstants->vz_fin_earth);

    // Get the initial position and velocity of the planet from startTime away from impact date.
    planet = planetInitial_incremental(0, startTime, planet, cConstants);
    

    // Setting the first element of planetCon to be equal to the planet conditions at startTime away from impact
    planetCon[0] = planet;

    // Shows progress of planet position calculations before the optimization in cuda can occur.
    std::cout << "Calculating planet positions for the trip time range" << std::endl;
    std::cout << "          10 20 30 40 50 60 70 80 90 100" << std::endl;
    std::cout << "progress:[";

    // Iterate backwards until total number of data points acquired
    for (int i = 1; i < tolData; i++) { 
        // Calculates the planet's condition at each point (time) from the previously calculated point.
        planetCon[i] = planetInitial_incremental(calc_time(i)-timeRes, calc_time(i), planet, cConstants); // Obtaining conditions of the planet

        // Filling progress bar, 30 value derived from the fact that there is 30 character spaces to fill in progress bar
        if ((i % (tolData/30)) == 0) {
            std::cout << ">";
        }

        // Sets planet equal to the conditions calculated for a given time to be used as reference point in next usage of planetInitial_incremental
        planet = planetCon[i];
    }
    // Closing progress bar and adding a couple of empty line spaces
    std::cout << "]\n\n";
}

// Returns conditions of the planet for a given time input
// Input: currentTime - time offset from impact backwards in time (larger value refers further back) in units of seconds
//          in other words, currentTimes how many seconds it will take for the spacecraft to reach the asteroid
// Output: Returns an element that is to the planet's position/velocity at currentTime away from impact using interpolate
//         as currentTime very likely does not directly corelate to an explicit derived element in planetCon
elements<double> PlanetInfo::getCondition(const double & currentTime) {
// __host__ __device__ elements<double> PlanetInfo::getCondition(const double & currentTime) {
    // Defining variables to use
    elements<double> lower;
    elements<double> upper;
    elements<double> result;
    
    int index;
    // Setting index equal to the nearest index of data based on time
    index = calcIndex(currentTime); 
    
    if (index < calcIndex(startTime)) {
        std::cout << "Planet condition index being accessed is less than startTime! Returning nearest index\n";
        index = calcIndex(startTime);
    }
    else if (index > calcIndex(endTime) - 1) {
        std::cout << "Planet condition index being accessed is greater than to endTime-1! Returning nearest valid index\n";
        index = calcIndex(endTime) - 1; // is (index-1) because the index variable is the lower of two, the keeps from going out of bounds with the upper index value 
    }

    // The lower index weight is equal to the index
    lower = planetCon[index];
    // The upper index weight is equal to the index plus 1
    upper = planetCon[index + 1];

    // Weights are used for interpolation
    double lowerWeight = 1 - ((currentTime-calc_time(index)) / timeRes);
    double upperWeight = 1 - ((calc_time(index+1)-currentTime) / timeRes);

    // Derive the element of the planet at this time using interpolate
    result = interpolate(lower,upper,lowerWeight,upperWeight);

    return result;
}

// Takes in a time and outputs a corresponding index (location of data).
int PlanetInfo::calcIndex(const double & currentTime) {
   return static_cast<int>((currentTime-startTime)/timeRes);
}

// Takes in an index and outputs the time corresponded to that index.
double PlanetInfo::calc_time(const int & currentIndex) {
    return startTime + (currentIndex*timeRes);
}

// Returns the tolData
int PlanetInfo::getTolData() {
    return tolData;
}

//make method to return pointer to array (planetCon)
elements<double>* PlanetInfo::getAllPositions(){
    return planetCon;
}

// Fills planet with the appropriate elements
void PlanetInfo::getPlanetPosition(const cudaConstants* cConstants, elements<double> & planet, const int & planetStatus){
    if(planetStatus == EARTH){
        planet = elements<double>(cConstants->r_fin_earth, cConstants->theta_fin_earth, cConstants->z_fin_earth, cConstants->vr_fin_earth, cConstants->vtheta_fin_earth, cConstants->vz_fin_earth);
    }
    if(planetStatus == MARS){
        planet = elements<double>(cConstants->r_fin_mars, cConstants->theta_fin_mars, cConstants->z_fin_mars, cConstants->vr_fin_mars, cConstants->vtheta_fin_mars, cConstants->vz_fin_mars);
    }
    
}



// Takes in the lower index, upper index, and their weights in order to calculate an approximate planet element for a time between two index.
// Input: lower & upper - elements that the interpolate is between of
//        lowerWeight & upperWeight - weight values to adjust the interpolation in cases where the time is closer to one index than another (not midpoint)
elements<double> PlanetInfo::interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight) {
    return (lower*lowerWeight)+(upper*upperWeight);
}

// Destructor, deallocates planetCon data
PlanetInfo::~PlanetInfo() {
    delete [] planetCon;
}

// Determines the next step back from a given element using rk4Reverse, used in the constructor of PlanetInfo
// Input for planetInitial_incremental:
//      timeInitial: time (in seconds) of the impact date
//      tripTime: optimized time period of the overall trip
//      planet: planet's position and velocity on the arrival date (passed in from planetInfo.cpp)
elements<double> planetInitial_incremental(double timeInitial, double tripTime, const elements<double> & planet, const cudaConstants * cConstants) {
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be finer than the time resolution.
  deltaT = (tripTime - timeInitial)/cConstants->min_numsteps; 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the planet's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial, tripTime, planet, deltaT, yp, cConstants->rk_tol, cConstants);
 
  return yp;
}

int getPlanetSize(const cudaConstants * cConstants){
    double startTime = 0.0; // Starting time (s), chronologically this is closest to impact time (0 would be exactly impact date)
    double endTime = cConstants->triptime_max;   // Ending time (s), chronologically this is earliest time away from impact date
    double timeRes = cConstants->timeRes;        // Time resolution for storing data points (s)
    int tolData = ((endTime-startTime)/timeRes) + 1; // Total Number of Data points in planetCon based on duration in seconds divided by resolution, plus one for the last 'section'
    return tolData*6;//size of array * 6 for elements
}


//interpolates the position of the planet currentTime seconds BEFORE the spacecraft reaches the asteroid
__host__ __device__ elements<double> getConditionDev(const double & currentTime, const cudaConstants * cConstants, const elements<double>* planetConditions) {
// __host__ __device__ elements<double> PlanetInfo::getCondition(const double & currentTime) {
    // Defining variables to use
    elements<double> lower;
    elements<double> upper;
    elements<double> result;
    double startTime = 0.0;
    double endTime = cConstants->triptime_max;
    double timeRes = cConstants->timeRes;
    
    int index, indexMin, indexMax;
    // Setting index equal to the nearest index of data based on time
    index = static_cast<int>((currentTime-startTime)/timeRes); 
    indexMin = 0;
    indexMax = static_cast<int>((endTime-startTime)/timeRes); 
    
    if (index < indexMin) {
        //std::cout << "Planet condition index being accessed is less than startTime! Returning nearest index\n";
        index = indexMin;
    }
    else if (index > indexMax- 1) {
        //std::cout << "Planet condition index being accessed is greater than to endTime-1! Returning nearest valid index\n";
        index = indexMax - 1; // is (index-1) because the index variable is the lower of two, the keeps from going out of bounds with the upper index value 
    }

    // The lower index weight is equal to the index
    lower = planetConditions[index];
    // The upper index weight is equal to the index plus 1
    upper = planetConditions[index + 1];

    // Weights are used for interpolation
    double lowerWeight = 1 - ((currentTime-(startTime + (index*timeRes))) / timeRes);
    double upperWeight = 1 - ((startTime + (((1+index)*timeRes)-currentTime)) / timeRes);

    // Derive the element of the planet at this time using interpolate
    result = (lower*lowerWeight)+(upper*upperWeight);

    return result;
}

