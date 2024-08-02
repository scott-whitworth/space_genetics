//Last Edit: August 2022
//Tasks Completed: 
	//Added gravitational assist from Mars
	//changes are in calcRate for vr, vtheta, vz
	
#include <math.h> // used for sine, cosine, and pow functions

template <class T>  __host__ __device__ elements<T> calc_k(const T & h, const elements<T>  & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal, const elements<T> & mars, const T & marsCraftDist) {
	return elements<T>( h*calcRate_r(y), 
						h*calcRate_theta(y), 
						h*calcRate_z(y), 
						h*calcRate_vr(y,coeff,accel,curTime, timeFinal, mars, marsCraftDist), 
						h*calcRate_vtheta(y,coeff,accel,curTime, timeFinal, mars, marsCraftDist),
						h*calcRate_vz(y,coeff,accel,curTime, timeFinal, mars, marsCraftDist));
}

template <class T> __host__ __device__ elements<T> calc_kPlanet(const T & h, const elements<T>  & y, const T & curTime, const T & timeFinal) {
	return elements<T>( h*calcRate_r(y), 
						h*calcRate_theta(y), 
						h*calcRate_z(y), 
						h*calcRate_vrPlanet(y), 
						h*calcRate_vthetaPlanet(y),
						h*calcRate_vzPlanet(y));
}

template <class T>  __host__ __device__ T calcRate_r(const elements<T> & y) {
	return y.vr;
}

template <class T>  __host__ __device__ T calcRate_theta(const elements<T> & y) {
	return y.vtheta / y.r;
}

template <class T>  __host__ __device__ T calcRate_z(const elements<T> & y) {
	return y.vz;
}

template <class T> __host__ __device__ T calcRate_vr(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal, const elements<T> & mars, const T & marsCraftDist) {
	return ((-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + 
		   ((-constG * massMars * (y.r - mars.r))/ (pow(marsCraftDist,3))) +
		   (pow(y.vtheta,2) / y.r) + 
		   (accel*cos(calc_tau(coeff,curTime, timeFinal))*sin(calc_gamma(coeff,curTime, timeFinal))));
}

template <class T> __host__ __device__ T calcRate_vtheta(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal, const elements<T> & mars, const T & marsCraftDist) {
	return  (static_cast<double>(0.0) + 
			((-constG * massMars *fmod(fmod(y.theta, 2*M_PI) - fmod(mars.theta, 2*M_PI), 2*M_PI) *((y.r+mars.r)/2)) / (pow(marsCraftDist,3))) +
			(-y.vr*y.vtheta / y.r) +
			(accel*cos(calc_tau(coeff,curTime, timeFinal))*cos(calc_gamma(coeff,curTime, timeFinal))));
}

template <class T> __host__ __device__ T calcRate_vz(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal, const elements<T> & mars, const T & marsCraftDist) {
	return	((-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), 1.5)) + 
			((-constG * massMars * (y.z - mars.z))  / (pow(marsCraftDist,3))) + 
		    static_cast<double>(0.0) +
		    (accel*sin(calc_tau(coeff,curTime, timeFinal))));
}
template <class T> __host__ __device__ T calcRate_vrPlanet(const elements<T> & y) {
	return ((-constG * massSun * y.r) / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + 
			(pow(y.vtheta,2) / y.r);
}

template <class T> __host__ __device__ T calcRate_vthetaPlanet(const elements<T> & y) {
	return (0.0) + 
			(-y.vr*y.vtheta / y.r);
}

template <class T> __host__ __device__ T calcRate_vzPlanet(const elements<T> & y) {
	return ((-constG * massSun * y.z) / pow(pow(y.r, 2) + pow(y.z, 2), 1.5)) + 
			(0.0);
}