#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>

//Non-configured constants

#define _USE_MATH_DEFINES // for use of M_PI
#define SECONDS_IN_YEAR (365.25*24*3600) // Used for converting units from years to seconds or vice-versa 31557600

//This values is often used to denote that we are not optimizing for a certain objective
#define NOT_APPLICABLE  -1

// Planetary properties and constants
#define earthRadius (1.49598261e11/AU) // radial distance of Earth from Sun (au)
#define earthPeriod SECONDS_IN_YEAR //Orbital period of the Earth around the Sun (s)
#define earthMass 5.97237e24 // mass of the earth (kg)
#define lunarMass 7.34767e22 // mass of the earth's moon (kg)
#define ESOI (earthRadius*pow(((earthMass+lunarMass)/massSun),0.4)) // escape sphere of influence (au)
#define AU 1.49597870691e11 // used to convert meters to astronomical units (m) 
#define constG 1.99349603314131e-44 // gravitational constant- used to calculate the gravitational force (AU^3/(s^2 * kg)) 
#define massSun 1.988500e30 // mass of the sun (kg)
#define massMars 6.41691e23 // mass of Mars (kg) -> source: https://ssd.jpl.nasa.gov/planets/phys_par.html#refs
#define marsPeriod 1.8808476*SECONDS_IN_YEAR //Orbital period of Mars around the Sun (s) -> source: https://ssd.jpl.nasa.gov/planets/phys_par.html#refs
#define marsRadius (3.38892E+06/AU) //radius of mars (au)
#define marsSunRadius (2.27940000000E+11/AU) //Distance of mars from the sun
#define MSOI (marsSunRadius*pow(((massMars)/massSun),0.4)) // escape sphere of influence (au)

// Starting location and sizes in the optimization array for navigation to access specific values
// Note: array sizes must be odd, coinciding with their use in computing Fourier series
// Reason why these 3 values are not in config is due to need for static expected memory with the GPU,
//        as well as current structural expectations for these values in compiling
//       (coefficients struct for example)
// WARNING: As of 7/6/2020, the file used in non-random start (optimizedVector.bin) is based on previous definition of values and so if these are changed from old sizes then it may cause issues for it
#define GAMMA_ARRAY_SIZE 7  // Number of coefficients for gamma
#define   TAU_ARRAY_SIZE 3  // Number of coefficients for tau
#define COAST_ARRAY_SIZE 5  // Number of coefficients for coasting

// Offset values and total length are based on the array sizes defined above
#define GAMMA_OFFSET     0                            // In the array of variables, GAMMA is at start index
#define TAU_OFFSET      (GAMMA_ARRAY_SIZE)            // TAU follows after GAMMA values
#define ALPHA_OFFSET    (TAU_OFFSET + TAU_ARRAY_SIZE) // ALPHA follows TAU
#define BETA_OFFSET     ( ALPHA_OFFSET + 1)           // BETA follows ALPHA
#define ZETA_OFFSET     (  BETA_OFFSET + 1)           // ZETA follows BETA
#define TRIPTIME_OFFSET (  ZETA_OFFSET + 1)           // TRIPTIME follows ZETA
#define COAST_OFFSET    (TRIPTIME_OFFSET + 1)         // COAST follows TRIPTIME

// OPTIM_VARS = Number of variables (array sizes plus 4 for alpha, beta, zeta, and triptime)
#define OPTIM_VARS (GAMMA_ARRAY_SIZE + TAU_ARRAY_SIZE + COAST_ARRAY_SIZE + 4)

//Constants used accessing info stored in the mission objectives vector for cConstants
#define MISSION_PARAMETER_OFFSET 0 //Where in the array the target parameter is stored
#define MISSION_OBJECTIVE_OFFSET 1 //Where in the array that the objective is stored
#define MISSION_GOAL_OFFSET 2      //Where in the array the target value is stored

//Enumerations for storing mission objectives within cConstants
//First enumeration is for identifying parameters
//  Will be used to grab the correct parameter from individuals
//TODO: as of August 2022, this is never used, recommend deleting this
enum MISSION_PARAMETERS {
    POS_DIFF = 0,       //The posDiff of an individual
    SPEED_DIFF = 1,     //The speedDiff of an individual
    FUEL_BURNED = 2,    //The amount of fuel burned by an individual
    TRIP_TIME = 3,      //The time an individual takes to complete its mission
};

//This is used in both Child and Adult
//enumeration to make error status easier to keep track of, as opposed to using hard-coded numbers 
//could be removed in the future, included for now to check if the default constructor is being used
//TODO: as of August 2022, this is currently not used outside of constructors, recommend deleting this
enum STATUS {
    DEFAULT_CHILD = 0, //child that is created through the default constructor, not ready to be made an adult
    FUNCTIONAL_CHILD = 1, //child that is created in main constructor and given startParams
    FUNCTIONAL_ADULT = 2, //once the child has gone through callRK and become an adult, we can change its status to FUNCTIONAL_ADULT
};

//status values used for child and adult that tell us if it is a nan and what kind of nan error it is
//    Addresses complicated communication from the GPU process that previously was paying attention only to NANs
//    Now the GPU loop can assign a status to more accurately describe why the individual is not VALID
enum ERROR_STATUS{
    VALID = 0,   //not a nan, no problems with any of the parameter values
    SUN_ERROR = 1,   //flew too close to the sun, the posDiff and speedDiff are set to bad values during callRK
    MARS_ERROR = 2, //the spacecraft gets too close to Mars
    BAD_ASSIST = 3, //The spacecraft had a net-negative orbital assist
    SIMNUM_ERROR = 4, //The individual needed more simulation cycles than the specified max
    DUPLICATE = 5,  //an adult has been found to be a duplicate of another, original adult
    NOT_RUN = 6, //child has not been run through callRK (this is the default when created)
    NAN_ERROR = 7, //there was an issue during the simulation of the child, causing its final pos to be nan
    OTHER_ERROR = 8,   //any nans not caught during callRk are set to this error status in optimization
};

//Enum will track the simulation status for each child
//  Each child will initially have a sim status of INITIAL_RUN
//  If a child enters an sphere of influence during it's run, it will have the status INSIDE_MSOI
//  If then the child exits the sphere of influence, it will be assigned OUTSIDE_SOI
//  This is used in the Runge-Kutta simulations to 
enum SIM_STATUS{
    INITIAL_SIM = 0,
    INSIDE_SOI = 1,
    OUTSIDE_SOI = 2,
    COMPLETED_SIM = 3,
};

enum PLANET{
    EARTH = 0,
    MARS = 1,
};

#define MAX_DISTANCE 1  //max distance for giveDistance, assigned to adults on the boundaries (AKA best and worst adults) - used to be 1.0e+12


#endif
