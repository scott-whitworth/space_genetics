// Didymos-Optimization_Project:

template <class T>
thruster<T>::thruster(const cudaConstants* gConfig) {
    // If no thruster, set values to 0
    if (gConfig->thruster_type == THRUST_TYPE::NO_THRUST) {
        type = THRUST_TYPE::NO_THRUST;
        P0 = 0;
    }

    // setting values (defined in header) for when type 1 is called (NEXT)
    else if (gConfig->thruster_type == THRUST_TYPE::NEXT_C) {
        type = THRUST_TYPE::NEXT_C;
        P0 = NEXTP0;
    }
    coastThreshold = gConfig->coast_threshold;
}

template <class T> T thruster<T>::calc_eff(const T & Pin) {
    // Data interpolation for thruster type NEXT
    if (type == THRUST_TYPE::NEXT_C) {
        return  -1.328086e-23*pow(Pin,6) + 6.207694e-19*pow(Pin,5) - 9.991813e-15*pow(Pin,4) +  7.701266e-11*pow(Pin,3) - 3.136031e-07*pow(Pin,2) +  6.805225e-04*Pin;       // Polynomial fit
    }
    else return 0;
}

template <class T> T thruster<T>::calc_m_Dot(const T & Pin) {
    if (type == THRUST_TYPE::NEXT_C) {
        if (Pin < 640) {
            return 0;
        }
        else if (Pin < 2550) {
            return 1.99E-06;
        }
        else if (Pin < 4500) {
            return 4.44E-06;
        }
        else {
            return NEXTm_Dot0;
        }
    }
    else return 0;
}


// template <class T> std::ostream & operator<<(std::ostream & Str, const thruster<T> & e) {
//     Str << std::fixed;
//     Str << std::setprecision(16); // number of decimals output into text file
//     Str << e.m_Dot << "\t" << e.P0 << "\n";
//     return Str;
// }

template <class T> __host__ __device__ T calc_accel(const T & radius, const T & z, thruster<T> & thrusterType, T & massExpelled, const T & deltaT, const bool & thrusting, const T & wetMass, const cudaConstants* cConstants) {

    if (cConstants->thruster_type == thruster<double>::NO_THRUST) {
        return 0;
    }

    // Thrusting is evaluated in calcFourier.cpp within calc_coast().
    // When thrusting is equal to zero, calc_accel() will not be evaluated.
    if (!thrusting) {
        return 0;
    }

    // If all of the fuel has been expelled, then no more thrust can be applied
    if (massExpelled >= cConstants->fuel_mass) {
        return 0;
    }

    // Defining variables for calc_accel().
    T Pin; // Power input
    T Pthrust; // Thrust power
    T thrust;

    // Power going into the spacecraft as a function of the radius of the spacecraft from the sun (r is non-dimensionalized by dividing by 1 AU).
    T separation = sqrt( pow(radius, 2) + pow(z, 2) );
    Pin = thrusterType.P0/pow(separation, 2); 

    //If the spacecraft is closer to the sun than the earth, the power in can not be greater than the power experimentally measured on earth.
    //This creates a "sphere" around the sun to ensure the power does not exceed the tested limit.
    if (separation <= 1) {
        Pin = thrusterType.P0; // It is devided by 1 astronomical unit to normalize it P0/(1 AU)
    }

    // The thrust power of the spacecraft is dependent upon the efficiency (calculated in thruster.cpp) and the power (in).
    Pthrust = thrusterType.calc_eff(Pin)*Pin; 

    // Thrust is calculated by power (thrust) and mDot.
    thrust = sqrt(2 * Pthrust * thrusterType.calc_m_Dot(Pin)); 

    // Calculates the amount of fuel used throughout the duration of the trip.
    massExpelled += thrusterType.calc_m_Dot(Pin) * deltaT;
    
    // the current mass of the spacecraft is equal to the fuel used minus the wetMass of the spacecraft
    // Acceleration of the spacecraft due to thrusting calculated by thrust divided by the mass of the spacecraft.
    // AU converts the acceleration from m/s^2 to au/s^2.
    return thrust/(AU*(wetMass - massExpelled));
}
