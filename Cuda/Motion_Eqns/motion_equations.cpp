//Didymos-Optimization_Project:
//Last Editor: Mateo and Lauren
//Tasks Completed: 
	//No recent changes
	
#include <math.h> // used for sine, cosine, and pow functions

template <class T>  __host__ __device__ elements<T> calc_k(const T & h, const elements<T>  & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal) {
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
						h*calcRate_vr(y,coeff,accel,curTime, timeFinal), h*calcRate_vtheta(y,coeff,accel,curTime, timeFinal),
						h*calcRate_vz(y,coeff,accel,curTime, timeFinal));
}

template <class T> __host__ __device__ elements<T> calc_kEarth(const T & h, const elements<T>  & y, const T & curTime, const T & timeFinal) {
	return elements<T>( h*calcRate_r(y), h*calcRate_theta(y), h*calcRate_z(y), 
						h*calcRate_vrEarth(y), h*calcRate_vthetaEarth(y),
						h*calcRate_vzEarth(y));
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

template <class T> __host__ __device__ T calcRate_vr(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal) {
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + (pow(y.vtheta,2) / y.r) +
		(accel*cos(calc_tau(coeff,curTime, timeFinal))*sin(calc_gamma(coeff,curTime, timeFinal)));
}

template <class T> __host__ __device__ T calcRate_vtheta(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal) {
	return -y.vr*y.vtheta / y.r + accel*cos(calc_tau(coeff,curTime, timeFinal))*cos(calc_gamma(coeff,curTime, timeFinal));
}

template <class T> __host__ __device__ T calcRate_vz(const elements<T> & y, coefficients<T> & coeff, const T & accel, const T & curTime, const T & timeFinal) {
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), 1.5)) + accel*sin(calc_tau(coeff,curTime, timeFinal));
}

template <class T> __host__ __device__ T calcRate_vrEarth(const elements<T> & y) {
	return (-constG * massSun * y.r / (pow(pow(y.r, 2) + pow(y.z, 2), 1.5))) + (pow(y.vtheta,2) / y.r);
}

template <class T> __host__ __device__ T calcRate_vthetaEarth(const elements<T> & y) {
	return -y.vr*y.vtheta / y.r;
}

template <class T> __host__ __device__ T calcRate_vzEarth(const elements<T> & y) {
	return (-constG * massSun * y.z / pow(pow(y.r, 2) + pow(y.z, 2), 1.5));
}