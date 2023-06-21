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

    // setting values (defined in header) for when type 1 is called (SPT_140)
    else if (gConfig->thruster_type == THRUST_TYPE::SPT_140) {
        type = THRUST_TYPE::SPT_140;
        P0 = SPTP0;     
    } 

    // setting values (defined in header) for when type 1 is called (AEPS)
    else if (gConfig->thruster_type == THRUST_TYPE::AEPS) {
        type = THRUST_TYPE::AEPS;
        P0 = AEPSP0;     
    }                                                 
    coastThreshold = gConfig->coast_threshold;
}

template <class T> T thruster<T>::calc_eff(const T & Pin) {
    // Data interpolation for thruster type NEXT
    if (type == THRUST_TYPE::NEXT_C) {
        return  -1.328086E-23*pow(Pin,6) + 6.207694E-19*pow(Pin,5) - 9.991813E-15*pow(Pin,4) +  7.701266E-11*pow(Pin,3) - 3.136031E-07*pow(Pin,2) +  6.805225E-04*Pin;       // Polynomial fit
    }
    else if (type == THRUST_TYPE::SPT_140) {
        return  -7.141057E-23*pow(Pin,6) + 1.703637E-18*pow(Pin,5) - 1.775691E-14*pow(Pin,4) + 1.035756E-10*pow(Pin,3) - 3.525689E-07*pow(Pin,2) + 6.628786E-04*Pin;   // Polynomial fit
    }
    else if (type == THRUST_TYPE::AEPS) {
        return   -1.869804E-24*pow(Pin,6) + 1.025302E-19*pow(Pin,5) - 2.221817E-15*pow(Pin,4) + 2.423269E-11*pow(Pin,3) - 1.398300E-07*pow(Pin,2) + 4.131622E-04*Pin;   // Polynomial fit
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
    else if (type == THRUST_TYPE::SPT_140) {
        if (Pin < 500){
            return 0;
        }
        else if (Pin < 6000) {
            return 3.044520E-09*Pin;
        }
        else {
            return SPT140m_Dot0;
        }
    }
    else if (type == THRUST_TYPE::AEPS) {
        if (Pin < 3120) {
            return 0;
        }
        else if (Pin < 6700) {
            return 3.248871E-09*Pin;
        }
        else  {
            return AEPSm_Dot0;
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
    T craftSunDist = sqrt( pow(radius, 2) + pow(z, 2) );
    Pin = thrusterType.P0/pow(craftSunDist, 2); 

    //If the spacecraft is closer to the sun than the earth, the power in can not be greater than the power experimentally measured on earth.
    //This creates a "sphere" around the sun to ensure the power does not exceed the tested limit.
    if (craftSunDist <= 1) {
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
