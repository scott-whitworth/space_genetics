//Didymos-Optimization_Project:
#include <iomanip> // setprecision(int)

// constructors

// sets starting values as given
template <class T>
elements<T>::elements(T r0, T theta0, T z0, T vr0, T vtheta0, T vz0) {
    r = r0;
    theta = theta0;
    z = z0;
    vr = vr0;
    vtheta = vtheta0;
    vz = vz0;
}

// sets values to default of zero
template <class T>
elements<T>::elements() {
    r = static_cast<T>(0);
    theta = static_cast<T>(0);
    z = static_cast<T>(0);
    vr = static_cast<T>(0);
    vtheta = static_cast<T>(0);
    vz = static_cast<T>(0);
}


//overload operators to do math on all the elements in the struct seperately
//Treating each element as a matrix operation
template <class T> elements<T> elements<T>::operator+(const elements<T> & e) const {
    return elements<T>(this->r + e.r, this->theta + e.theta, this->z + e.z, this->vr + e.vr, this->vtheta + e.vtheta, this->vz + e.vz);
}

template <class T> elements<T> elements<T>::operator-(const elements<T> & e) const {
    return elements<T>(this->r - e.r, this->theta - e.theta, this->z - e.z, this->vr - e.vr, this->vtheta - e.vtheta, this->vz - e.vz);
}

// template <class T> elements<T> elements<T>::operator*(const elements<T> & e) const {
//     return elements<T>(this->r * e.r, this->theta * e.theta, this->z * e.z, this->vr * e.vr, this->vtheta * e.vtheta, this->vz * e.vz);
// }

// template <class T> elements<T> elements<T>::operator/(const elements<T> & e) const {
//     return elements<T>(this->r / e.r, this->theta / e.theta, this->z / e.z, this->vr / e.vr, this->vtheta / e.vtheta, this->vz / e.vz);
// }

//constructor which takes in an scalar
template <class T> elements<T> elements<T>::operator*(const T & i) const {
    return elements<T>(this->r * i,  this->theta * i, this->z * i, this->vr * i, this->vtheta * i, this->vz * i);
}

template <class T> elements<T> elements<T>::operator/(const T & i) const {
    return elements<T>( this->r / i, this->theta / i, this->z / i, this->vr / i, this->vtheta / i, this->vz / i);
}

// Comparison method to compare this elements to another, returning true (equivalent) within a threshold
// input: other - another elements object to compare to
// output: returns the percentage on how different the two elements are
template <class T> double elements<T>::compare(const elements<T> & other) {
    //Check values
    double percentDiff = 0;
    percentDiff += abs(100*abs(this->r - other.r)/(0.5*(this->r + other.r)));
    percentDiff += abs(100*abs(this->theta - other.theta)/(0.5*(this->theta + other.theta)));
    percentDiff += abs(100*abs(this->z - other.z)/(0.5*(this->z + other.z)));
    percentDiff += abs(100*abs(this->vr - other.vr)/(0.5*(this->vr + other.vr)));
    percentDiff += abs(100*abs(this->vtheta - other.vtheta)/(0.5*(this->vtheta + other.vtheta)));
    percentDiff += abs(100*abs(this->vz - other.vz)/(0.5*(this->vz + other.vz)));

    return percentDiff;
}

template <class T> std::ostream & operator<<(std::ostream & Str, const elements<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into .csv file
    Str << std::scientific;
    Str << e.r << "," << e.theta << "," << e.z << "," << e.vr << "," << e.vtheta << "," << e.vz << "\n";
    return Str;
}
