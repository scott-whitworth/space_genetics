template <class T> rkParameters<T>::rkParameters(T tripTime0,T r0, T theta0, T z0, T vr0, T vtheta0, T vz0, // elements<T>
                                                 T *gamma0, T *tau0, T *coast0) {  // coefficients<T>
    tripTime = tripTime0;

    y0.r = r0;
    y0.theta = theta0;
    y0.z = z0;
    y0.vr = vr0;
    y0.vtheta = vtheta0;
    y0.vz = vz0;

    for (int i = 0; i < coeff.gammaSize; i++) {
        coeff.gamma[i] = gamma0[i];
    }
    for (int i = 0; i < coeff.tauSize; i++) {
        coeff.tau[i] = tau0[i];
    }
    for (int i = 0; i < coeff.coastSize; i++) {
        coeff.coast[i] = coast0[i];
    }
}

template <class T> rkParameters<T>::rkParameters(T tripTime0,  T r0, T theta0, T z0, T vr0, T vtheta0, T vz0,T *gamma0, T *tau0, T *coast0, T alpha0, T beta0, T zeta0) {
    
    tripTime = tripTime0;

    y0.r = r0;
    y0.theta = theta0;
    y0.z = z0;
    y0.vr = vr0;
    y0.vtheta = vtheta0;
    y0.vz = vz0;

    for(int i = 0; i < coeff.gammaSize; i++){
        coeff.gamma[i] = gamma0[i];
    }
    for(int i = 0; i < coeff.tauSize; i++){
        coeff.tau[i] = tau0[i];
    }
    for(int i = 0; i < coeff.coastSize; i++){
        coeff.coast[i] = coast0[i];
    }

    alpha = alpha0;
    beta = beta0;
    zeta = zeta0;
} 

template <class T> rkParameters<T>::rkParameters(T tripTime0, elements<T> initialCondition, coefficients<T> coeff0) {
    tripTime = tripTime0;

    y0 = initialCondition;
    coeff = coeff0;
}

template <class T> rkParameters<T>::rkParameters(T tripTime0, T alpha0, T beta0, T zeta0, coefficients<T> coeff0) {
    tripTime = tripTime0;
    alpha = alpha0;
    beta = beta0;
    zeta = zeta0;
    coeff = coeff0;

    y0.r = 0;
    y0.theta = 0;
    y0.z = 0;
    y0.vr = 0;
    y0.vtheta = 0;
    y0.vz = 0;
} 

template <class T> rkParameters<T>::rkParameters() {
    tripTime = 0;

    y0.r = 0;
    y0.theta = 0;
    y0.z = 0;
    y0.vr = 0;
    y0.vtheta = 0;
    y0.vz = 0;

    alpha = 0;
    beta = 0;
    zeta = 0;

    //causes error: "expression must be a modifiable lvalue"
    //coeff.gamma = NULL;
    //coeff.tau = NULL;
    //coeff.coast = NULL;
}

// Method for creating random parameters
// input: rng - random number object generator to be used to generate random values
//        cConstants - to access the random range values
// output: an rkParameters object that contains randomized properties within a valid ranges
rkParameters<double> randomParameters(std::mt19937_64 & rng, const cudaConstants * cConstants) {
    double tripTime =  ( ((cConstants->triptime_max - cConstants->triptime_min) * (static_cast<double>(rng()) / rng.max())) + cConstants->triptime_min);

    double alpha = cConstants->alpha_random_start_range * 2 * ((static_cast<double>(rng()) / rng.max()) - 0.5);
    double  beta = cConstants-> beta_random_start_range * ((static_cast<double>(rng()) / rng.max())); // From 0 to start range, beta cannot be negative
    double  zeta = cConstants-> zeta_random_start_range * 2 * ((static_cast<double>(rng()) / rng.max()) - 0.5); 

    coefficients<double> testcoeff;
    for (int j = 0; j < testcoeff.gammaSize; j++) {
        testcoeff.gamma[j] = cConstants->gamma_random_start_range * 2 * ((static_cast<double>(rng()) / rng.max()) - 0.5);
    }
    for (int j = 0; j < testcoeff.tauSize; j++) {
        testcoeff.tau[j] = cConstants->tau_random_start_range * 2 * ((static_cast<double>(rng()) / rng.max()) - 0.5);
    }
    for (int j = 0; j < testcoeff.coastSize; j++) {
        testcoeff.coast[j] = cConstants->coast_random_start_range * 2 * ((static_cast<double>(rng()) / rng.max()) - 0.5);
    }

    //Initialize the output parameter
    return rkParameters<double>(tripTime, alpha, beta, zeta, testcoeff); 
}


template <class T> bool rkParameters<T>::compare(const rkParameters<T> & other, T comp_Thresh) {
    //First Check Coefficient element
    for (int i = 0; i < this->coeff.gammaSize; i++) {
        if (abs(this->coeff.gamma[i] - other.coeff.gamma[i]) > comp_Thresh) {
            return false;
        }
    }
    for (int i = 0; i < this->coeff.tauSize; i++) {
        if (abs(this->coeff.tau[i] - other.coeff.tau[i]) > comp_Thresh) {
            return false;
        }
    }
    for (int i = 0; i < this->coeff.coastSize; i++) {
        if (abs(this->coeff.coast[i] - other.coeff.coast[i]) > comp_Thresh) {
            return false;
        }
    }

    //Check Starting pos/vel
    if ( !this->y0.compare(other.y0,comp_Thresh) ) {
        return false;
    }

    //Other Conditions
    if (abs(this->tripTime - other.tripTime) > comp_Thresh) {
        return false;
    }
    /*
    if( abs(this->wetMass - other.wetMass) > comp_Thresh){
        return false;
    }
    */

    //If we have made it this far, everthing is good
    return true;
}

template <class T> void rkParameters<T>::parametersRK4Simple(T timeInitial, T stepSize, T absTol, elements<T> & y) {
    double accel = 0;
    rk4Simple(timeInitial, tripTime, y0, stepSize, y, absTol, coeff, accel, static_cast<double>(WET_MASS));
}

template <class T> std::ostream & operator<<(std::ostream & Str, const rkParameters<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    
    Str << "Coeff: \n" << e.coeff;
    Str << "Elements:\n\t" << e.y0;
    Str << "Final Time: " << e.tripTime << std::endl;
    Str << "Launch angles:\n\tAlpha: " << e.alpha << "\n\tBeta: " << e.beta << "\n\tZeta: " << e.zeta << std::endl;
    //" WetMass: " << e.wetMass << endl;
    return Str;
}