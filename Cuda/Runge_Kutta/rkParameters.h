#ifndef RKPARAMETERS_H
#define RKPARAMETERS_H

//structs
#include "../Motion_Eqns/elements.h"
#include "../Thrust_Files/coefficients.h"

//struct to hold all the values required for the runge-kutta functions
template <class T> struct rkParameters {
    
    //////////////////
    // Constructors //
    //////////////////

    // Updated constructor which sets all the components according to values taken in
   __host__ __device__ rkParameters<T>(T tripTime0, 
                   T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                   T *gamma0, T *tau0, T *coast0, T alpha0, T beta0, T zeta0); // coefficients<T>


    // Constructor which sets all the components according to values taken in
   __host__ __device__ rkParameters<T>(T tripTime0, 
                   T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                   T *gamma0, T *tau0, T *coast0); // coefficients<T>

    // Alternate Constructor
    // timeFinal0           - Total time of mission (s)
    // initialCondition     - Initial position of the spacecraft at start of calculation (escape position/velocity)
    // coeff0               - Coefficients for the thrust angles and  acceleration
    __host__ __device__ rkParameters<T>(T tripTime0,
                  elements<T> initialCondition, coefficients<T> coeff0);  

    // Updated alternate Constructor
    // timeFinal0           - Total time of mission (s)
    // coeff0               - Coefficients for the thrust angles and  acceleration
    // alpha,beta,zeta      - Launch angles at SOI 
    __host__ __device__ rkParameters<T>(T tripTime0, T alpha0, T beta0, T zeta0, coefficients<T> coeff0);  



    // constructor which sets everything to zero
    __host__ __device__ rkParameters<T>();
    
    /////////////
    // Members //
    /////////////

    // Dependent Variables:

    // Initial Position/Velocity elements
    // Contains r, theta, z, vr, vtheta, and vz
    elements<T> y0;

    // Optimized Variables:

    // Initial Optimization Coefficients
    // Contains arrays for Fourier series:
    //    gamma
    //    tau
    //    coasting
    // and a value for coast_threshold
    coefficients<T> coeff;

    // Final Time of simulation (s)
    T tripTime;

    // launch angles at SOI
    T alpha;    // position
    T beta;     // velocity in-plane
    T zeta;     // velocity out-of-plane

    /////////////////////
    // Utility Methods //
    /////////////////////
    // Comparison function (for sort())
    // Param: other - another rkParameter to be compared to
    // Param: comp_Thresh - comparison threshold
    // Returns true all elements of other are the same as *this, within the threshold comp_Thresh
    // Not exact, as there are a number of different magnitudes represented in parameters
    bool compare(const rkParameters<T> & other, T comp_Thresh);

    // Based on *this parameters, calculates the final position of the spacecraft using CPU based Methods
    // runs the CPU version of the RK calculation to find the final position of the spacecraft, using rkParameters for input
    // input the values which are not run-specific
    // Param: timeInitial - start time of simulation (s)
    // Param: stepSize - first time interval between data points (s)
    // Param: absTol - error tolerence for Runge-Kutta
    // Output: y will contain the final position of the simulation, y initial value will be overwritten
    void parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y);

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const rkParameters<T> & e); 
};

// Method for creating random parameters
// input: rng - random number object generator to be used to generate random values
//        cConstants - to access the random range values
// output: an rkParameters object that contains randomized properties within a valid ranges
rkParameters<double> randomParameters(std::mt19937_64 & rng, const cudaConstants * cConstants);

#include "rkParameters.cpp"
#endif