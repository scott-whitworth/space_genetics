#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>

#define _USE_MATH_DEFINES // for use of M_PI
#define SECONDS_IN_YEAR (365.25*24*3600) // Used for converting units from years to seconds or vice-versa 31557600

// Planetary properties and constants
#define earthRadius (1.49598261e11/AU) // radial distance of Earth from Sun (au)
#define earthPeriod SECONDS_IN_YEAR //Orbital period of the Earth around the Sun (s)
#define earthMass 5.97237e24 // mass of the earth (kg)
#define lunarMass 7.34767e22 // mass of the earth's moon (kg)
#define ESOI (earthRadius*pow(((earthMass+lunarMass)/massSun),0.4)) // escape sphere of influence (au)
#define AU 1.49597870691e11 // used to convert meters to astronomical units (m) 
#define constG 1.99349603314131e-44 // gravitational constant- used to calculate the gravitational force (AU^3/(s^2 * kg)) 
#define massSun 1.988500e30 // mass of the sun (kg)

//Definitions for mission type
#define Rendezvous 1
#define Impact 2

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


//status values used for child and adult that tell us if it is a nan and what kind of nan error it is
#define VALID       0   //not a nan, no problems with any of the parameter values
#define SUN_ERROR   1   //flew too close to the sun, the posDiff and speedDiff are set to bad values during callRK
#define OTHER_ERROR 2   //any nans not caught during callRk are set to this error status in optimization

//Error values for nans to be changed to, used in callRK and optimization when finding nans
#define BAD_POSDIFF          10 //a bad posDiff for either mission to have (10 AU)
#define BAD_RENDEV_SPEEDDIFF 10 //a bad speedDiff for a rendezvous mission to have (10 AU/s)
#define BAD_IMPACT_SPEEDDIFF 0  //a bad speedDiff for an impact mission to have (0 AU/s)


#endif
