% impactParams(): script to calculate position and velocity vector
% components in cylindrical coordinates at impact using Cartesian vector
% components retrieved from https://ssd.jpl.nasa.gov/sbdb.cgi 
% impact date: 09-30-2022 19:54:55 UTC (data available at 19:55:00)
clc,clear
day2sec = 3600*24;

% Asteroid position
x_a = 6.061727293707692E-01; y_a = 8.369531060353675E-01;
z_a = -2.015445997254413E-02;
[th_a, r_a] = cart2pol(x_a,y_a);
r_a = vpa(r_a), th_a = vpa(th_a), z_a = vpa(z_a)
% Asteroid velocity
vx_a = -1.472521774199771E-02/day2sec; vy_a = 1.321610726487941E-02/day2sec;
vz_a = 1.066836745240369E-03/day2sec;
vr_a = (x_a*vx_a+y_a*vy_a)/r_a; 
vth_a = (x_a*vy_a-y_a*vx_a)/r_a; 
vr_a = vpa(vr_a), vth_a = vpa(vth_a), vz_a = vpa(vz_a)

% Earth position
x_e = 9.919106263356792E-01; y_e = 1.363370226016378E-01;
z_e = -9.245804016889583E-06;
[th_e, r_e] = cart2pol(x_e,y_e);
r_e = vpa(r_e), th_e = vpa(th_e), z_e = vpa(z_e)
% Earth velocity
vx_e = -2.622869688402269E-03/day2sec; vy_e = 1.697947317059088E-02/day2sec;
vz_e = -7.867184752439903E-07/day2sec;
vr_e = (x_e*vx_e+y_e*vy_e)/r_e; 
vth_e = (x_e*vy_e-y_e*vx_e)/r_e; 
vr_e = vpa(vr_e), vth_e = vpa(vth_e), vz_e = vpa(vz_e)