
#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI

#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y_new, 
                               const T & absTol, coefficients<T> coeff, T *gamma,  T *tau, int & lastStep, T *accel_output, T *fuelSpent, const T & wetMass, const cudaConstants* cConstant) {

    thruster<T> thrust(cConstant);

    T curTime = timeInitial; // setting time equal to the start time
    int n = 0; // setting the initial iteration number equal to 0

    //mass of fuel expended (kg)
    //set to 0 initially
    T massFuelSpent = 0;

    // u - current position
    // error needs to be defined, but not being used
    elements<T> u, error;

    // Set the first element of the solution vector to the initial conditions
    u = y0;

    while (curTime <= timeFinal) { // iterate until time is equal to the stop time

        y_new[n] = u;

        times[n] = curTime;

        // Check the thruster type before performing calculations, at time 0
        if (cConstant->thruster_type == thruster<double>::NO_THRUST) {
            gamma[n] = tau[n] = accel_output[n] = fuelSpent[n] = 0;
        }
        else {
            // array of gamma for binary output
            gamma[n] = calc_gamma(coeff, curTime, timeFinal);
            // array of tau for binary output
            tau[n] = calc_tau(coeff,curTime, timeFinal);
            // array of acceleration for binary output
            accel_output[n] = calc_accel(u.r,u.z, thrust, massFuelSpent, stepSize, calc_coast(coeff, curTime, timeFinal, thrust), wetMass, cConstant);
            // array of fuel spent for binary output
            fuelSpent[n] = massFuelSpent;
        }

        if (curTime == timeFinal) {
            break;
        }
        
        //calculate new position
        rkCalc(curTime, timeFinal, stepSize, u, coeff, accel_output[n], error);

        //array of time output as t         
        curTime += stepSize;

        // Choosing a constant max number of steps for high precision final output
        stepSize = (timeFinal-timeInitial) / cConstant->cpu_numsteps;
        
        // shorten the last step to end exactly at time final
        if ( (curTime+stepSize) > timeFinal) {
            stepSize = (timeFinal-curTime);
        }

        n++;
    } //end of while 
    lastStep = n;
    
    // Test outputs to observe difference between rk4sys results with CUDA runge-kutta results
    std::cout << "rk4sys posDiff: " << sqrt(pow(cConstant->r_fin_ast - y_new[lastStep].r, 2) + pow(cConstant->r_fin_ast * cConstant->theta_fin_ast - y_new[lastStep].r * fmod(y_new[lastStep].theta, 2 * M_PI), 2) + pow(cConstant->z_fin_ast - y_new[lastStep].z, 2)) << std::endl;
    std::cout << "rk4sys velDiff: " << sqrt(pow(cConstant->vr_fin_ast - y_new[lastStep].vr, 2) + pow(cConstant->vtheta_fin_ast - y_new[lastStep].vtheta, 2) + pow(cConstant->vz_fin_ast - y_new[lastStep].vz, 2));
}

// ** Currently not used **
template <class T> void rk4Simple(const T & timeInitial, const T & timeFinal, const elements<T> & y0,
                                    T stepSize, elements<T> & y_new, const T & absTol, coefficients<T> coeff, T & accel, const T & wetMass, const cudaConstants * cConstants) {
    // Set the first element of the solution vector to the initial conditions of the spacecraft
    y_new = y0;

    thruster<T> thrust(cConstants);

    T curTime = timeInitial; // setting time equal to the start time

    //mass of fuel expended (kg)
    //set to 0 initially
    T massFuelSpent = 0;

    // Setting step size to be used until the last step
    stepSize = (timeFinal-timeInitial) / cConstants->cpu_numsteps;
    
    elements<T> error;
    bool coast;
    
    while (curTime < timeFinal) {  // iterate until time is equal to the stop time

        // Check the thruster type before performing calculations
        if (cConstants->thruster_type == thruster<double>::NO_THRUST) {
            coast = accel = 0;
        }
        else {
            // defining coast using calc_coast()
            coast = calc_coast(coeff, curTime, timeFinal, thrust);
            // defining acceleration using calc_accel()
            accel = calc_accel(y_new.r,y_new.z, thrust, massFuelSpent, stepSize, coast, wetMass);
        }

        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y_new, coeff, accel, error); 

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol, cConstants->doublePrecThresh);

        //The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize > (timeFinal-timeInitial) / cConstants->min_numsteps) {
            stepSize = (timeFinal-timeInitial) / cConstants->min_numsteps;
        }
        else if (stepSize < ((timeFinal-timeInitial) / cConstants->max_numsteps)) {
            stepSize = (timeFinal-timeInitial) / cConstants->max_numsteps;
        }

        // shorten the last step to end exactly at time final
        if ((curTime+stepSize) > timeFinal) {
            stepSize = (timeFinal-curTime);
        }

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft increases to 1000, so that path is not used for optimization.
        if (sqrt( pow(y_new.r,2) + pow(y_new.z,2) ) < 0.5) {
            y_new.r = 1000;
            return;
        }
        
    }  //end of while 
}

template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
                                   T stepSize, elements<T> & y_new, const T & absTol, const cudaConstants * cConstants) {
    // Set the first element of the solution vector to the conditions of earth on impact date (Oct. 5, 2022)
    y_new = y0;
    elements<T> error;
    T curTime = timeFinal; // setting time equal to the start time

    while( curTime > timeInitial) {  // iterates in reverse
        //calculate k values
        rkCalcEarth(curTime, timeFinal, stepSize, y_new, error);

        //array of time output as t         
        curTime -= stepSize;

        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol, cConstants->doublePrecThresh);

        //The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize > (timeFinal-timeInitial) / cConstants->min_numsteps) {
            stepSize = (timeFinal-timeInitial) / cConstants->min_numsteps;
        }
        else if (stepSize < ((timeFinal-timeInitial) / cConstants->max_numsteps)) {
            stepSize = (timeFinal-timeInitial) / cConstants->max_numsteps;
        }

        // shorten the last step to end exactly at time final
        if ( (curTime-stepSize) < timeInitial) {
            stepSize = curTime-timeInitial;
        }
    } //end of while
}

template <class T> __host__ __device__ void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, coefficients<T> & coeff, const T & accel, 
                                                    elements<T> & error) {

    // k variables for Runge-Kutta calculation of y_new
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    // Coefficients from MATLAB's implementation of ode45
    // Our calculation of k has the time step built into it (see motion_equations.cpp)
    k1 = calc_k(stepSize, y_new, coeff, accel, curTime, timeFinal); 
    k2 = calc_k(stepSize, y_new+k1*(static_cast <double> (1)/static_cast <double> (5)), coeff, accel, curTime+((static_cast <double> (1)/static_cast <double> (5))*stepSize), timeFinal); 
    k3 = calc_k(stepSize, y_new+k1*(static_cast <double> (3)/static_cast <double> (40))+k2*(static_cast <double> (9)/static_cast <double> (40)), coeff, accel, curTime+((static_cast <double> (3)/static_cast <double> (10))*stepSize), timeFinal);   
    k4 = calc_k(stepSize, y_new+k1*(static_cast <double> (44)/static_cast <double> (45))+k2*(static_cast <double> (-56)/static_cast <double> (15))+k3*(static_cast <double> (32)/static_cast <double> (9)), coeff, accel, curTime+((static_cast <double> (4)/static_cast <double> (5))*stepSize), timeFinal); 
    k5 = calc_k(stepSize, y_new+k1*(static_cast <double> (19372)/static_cast <double> (6561))+k2*(static_cast <double> (-25360)/static_cast <double> (2187))+k3*(static_cast <double> (64448)/static_cast <double> (6561))+k4*(static_cast <double> (-212)/static_cast <double> (729)), coeff, accel, curTime+((static_cast <double> (8)/static_cast <double> (9))*stepSize), timeFinal); 
    k6 = calc_k(stepSize, y_new+k1*(static_cast <double> (9017)/static_cast <double> (3168))+k2*(static_cast <double> (-355)/static_cast <double> (33))+k3*(static_cast <double> (46732)/static_cast <double> (5247))+k4*(static_cast <double> (49)/static_cast <double> (176))+k5*(static_cast <double> (-5103)/static_cast <double> (18656)), coeff, accel, curTime+stepSize, timeFinal);  
    k7 = calc_k(stepSize, y_new+k1*(static_cast <double> (35)/static_cast <double> (384))+k3*(static_cast <double> (500)/static_cast <double> (1113))+k4*(static_cast <double> (125)/static_cast <double> (192))+k5*(static_cast <double> (-2187)/static_cast <double> (6784))+k6*(static_cast <double> (11)/static_cast <double> (84)), coeff, accel, curTime+stepSize, timeFinal);  

    // New value
    y_new = y_new + k1*(static_cast <double> (35)/static_cast <double> (384)) + k3*(static_cast <double> (500)/static_cast <double> (1113)) + k4*(static_cast <double> (125)/static_cast <double> (192)) - k5*(static_cast <double> (2187)/static_cast <double> (6784)) + k6*(static_cast <double> (11)/static_cast <double> (84)) + k7*(static_cast <double> (0)/static_cast <double> (40));  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    //y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    //error = y_new-y_prev;

    // MATLAB code : ERROR between GPU and CPU
    //error = k1*(71)/(57600) + k3*(-71)/(16695) + k4*(71)/(1920)
    //- k5*(17253)/(339200) + k6*(22)/(525) + k7*(-1)/(40);

    // (outdated comment, but still important) Without k7 : no error between GPU and CPU (this has to do with the comment below)
    
    // Comonents of error are going to be really small. Need to make sure they are not too small to do anything with in calc_scalingFactor
    error =  ((k1*(static_cast <double> (71)/static_cast <double> (57600))) + (k3*(static_cast <double> (-71)/static_cast <double> (16695))) + (k4*(static_cast <double> (71)/static_cast <double> (1920)))  + (k5*(static_cast <double> (-17253)/static_cast <double> (339200))) + (k6*(static_cast <double> (22)/static_cast <double> (525)))) + (k7*(static_cast <double> (-1)/static_cast <double> (40)));
}

// The stepSize value that is inputted is assumed to be a positive value
template <class T> void rkCalcEarth(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error) {
    // Runge-Kutta algorithm    
    // k variables for Runge-Kutta calculation of y_new
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    
    stepSize *= -1; // Make this copy of stepSize negative as it goes backwards

    //calc_k multiplies all values by the stepSize internally.
    k1 = calc_kEarth(stepSize, y_new, curTime, timeFinal);        
    k2 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (1)/static_cast <double> (5))), curTime+((static_cast <double> (1)/static_cast <double> (5))*stepSize), timeFinal);   
    k3 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (3)/static_cast <double> (40)))+(k2*(static_cast <double> (9)/static_cast <double> (40))), curTime+((static_cast <double> (3)/static_cast <double> (10))*stepSize), timeFinal);   
    k4 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (44)/static_cast <double> (45)))+(k2*(static_cast <double> (-56)/static_cast <double> (15)))+(k3*(static_cast <double> (32)/static_cast <double> (9))), curTime+((static_cast <double> (4)/static_cast <double> (5))*stepSize), timeFinal);    
    k5 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (19372)/static_cast <double> (6561)))+(k2*(static_cast <double> (-25360)/static_cast <double> (2187)))+(k3*(static_cast <double> (64448)/static_cast <double> (6561)))+(k4*(static_cast <double> (-212)/static_cast <double> (729))), curTime+((static_cast <double> (8)/static_cast <double> (9))*stepSize), timeFinal);        
    k6 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (9017)/static_cast <double> (3168)))+(k2*(static_cast <double> (-355)/static_cast <double> (33)))+(k3*(static_cast <double> (46732)/static_cast <double> (5247)))+(k4*(static_cast <double> (49)/static_cast <double> (176)))+(k5*(static_cast <double> (-5103)/static_cast <double> (18656))), curTime+stepSize, timeFinal);        
    k7 = calc_kEarth(stepSize, y_new+(k1*(static_cast <double> (35)/static_cast <double> (384)))+(k3*(static_cast <double> (500)/static_cast <double> (1113)))+(k4*(static_cast <double> (125)/static_cast <double> (192)))+(k5*(static_cast <double> (-2187)/static_cast <double> (6784)))+(k6*(static_cast <double> (11)/static_cast <double> (84))), curTime+stepSize, timeFinal);  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    //v = y_new + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + (k1*(static_cast <double> (35)/static_cast <double> (384))) + (k3*(static_cast <double> (500)/static_cast <double> (1113))) + (k4*(static_cast <double> (125)/static_cast <double> (192))) + (k5*(static_cast <double> (-2187)/static_cast <double> (6784))) + (k6*(static_cast <double> (11)/static_cast <double> (84)));  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    //y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    //error = y_new-y_prev;

    error = (k1*(static_cast <double> (71)/static_cast <double> (57600))) + (k3*(static_cast <double> (-71)/static_cast <double> (16695))) + (k4*(static_cast <double> (71)/static_cast <double> (1920))) - (k5*(static_cast <double> (17253)/static_cast <double> (339200))) + (k6*(static_cast <double> (22)/static_cast <double> (525))) + (k7*(static_cast <double> (-1)/static_cast <double> (40)));    
}

template <class T> __host__ __device__ T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, const double precThresh) {
    // relative total error is the total error of all coponents of y which is used in scale.
    // scale is used to determine the next step size.
    T normTotError, scale;

    // relative error (unitless) 
    elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
    difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);

    if (!pmLimitCheck(pmError, precThresh)) {
        // pmError is too small!
        // Keep time step the same
        // Complicated rational:
        // If error is so small the error is at the double error limit
        //     - Could mean time step needs to be smaller, thus more presice
        //     - Or could mean that time step needs to be longer to get a reasonable error
        return 1.0; // Change step size by 0%
    }

    // elements<T> pmError(previous.r, previous.theta, previous.z, previous.vr,  previous.vtheta, previous.vz);

    // square root of sum of squares of the error from the 6 elements to determine the scale for the time step of the next iteration
    normTotError = sqrt(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2));
    scale = pow((absTol/normTotError),0.2);

    return scale;   
}

template <class T> __host__ __device__ bool pmLimitCheck(const elements<T> & pmError, const double precThresh){
    //It is possible this is a major resource drain. This might be faster to square everything and not use fabs (floating point abs)
    if( (fabs(pmError.r) < precThresh ) ||
        (fabs(pmError.theta) < precThresh ) ||
        (fabs(pmError.z) < precThresh ) ||
        (fabs(pmError.vr) < precThresh ) ||
        (fabs(pmError.vtheta) < precThresh ) ||
        (fabs(pmError.vz) < precThresh ) )
    {
        //Error it too small for precise calculation of step size
        return false;
    } else {
        //Error is large enough that accurate step size can be computed
        // All error values are at least within 14 orders of magnitued of their original guesses
        return true;
    }
}