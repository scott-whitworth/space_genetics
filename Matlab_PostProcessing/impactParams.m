% impactParams(): script to calculate position and velocity vector
% components in cylindrical coordinates at impact using Cartesian vector
% components retrieved from https://ssd.jpl.nasa.gov/sbdb.cgi 
% Arrival date: 01-26-2026 00:00:00 UTC (data available at 00:00:00)
clc,clear

params = xlsread('Psyche_impactParams.xlsx');
day2sec = 3600*24;

% Earth position
x_e = params(1); y_e = params(2);
z_e = params(3);
[th_e, r_e] = cart2pol(x_e,y_e);
r_e = vpa(r_e), th_e = vpa(th_e), z_e = vpa(z_e)
% Earth velocity
vx_e = params(4')/day2sec; vy_e = params(5)/day2sec;
vz_e = params(6)/day2sec;
vr_e = (x_e*vx_e+y_e*vy_e)/r_e; 
vth_e = (x_e*vy_e-y_e*vx_e)/r_e; 
vr_e = vpa(vr_e), vth_e = vpa(vth_e), vz_e = vpa(vz_e)

% Target position
x_t = params(7); y_t = params(8);
z_t = params(9);
[th_t, r_t] = cart2pol(x_t,y_t);
r_t = vpa(r_t), th_t = vpa(th_t), z_t = vpa(z_t)
% Target velocity
vx_t = params(10)/day2sec; vy_t = params(11)/day2sec;
vz_t = params(12)/day2sec;
vr_t = (x_t*vx_t+y_t*vy_t)/r_t; 
vth_t = (x_t*vy_t-y_t*vx_t)/r_t; 
vr_t = vpa(vr_t), vth_t = vpa(vth_t), vz_t = vpa(vz_t)

% Mars position
x_m = params(13); y_m = params(14);
z_m = params(15);
[th_m, r_m] = cart2pol(x_m,y_m);
r_m = vpa(r_m), th_m = vpa(th_m), z_m = vpa(z_m)
% Mars velocity
vx_m = params(16)/day2sec; vy_m = params(17)/day2sec;
vz_m = params(18)/day2sec;
vr_m = (x_m*vx_m+y_m*vy_m)/r_m; 
vth_m = (x_m*vy_m-y_m*vx_m)/r_m; 
vr_m = vpa(vr_m), vth_m = vpa(vth_m), vz_m = vpa(vz_m)