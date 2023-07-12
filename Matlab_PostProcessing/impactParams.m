% impactParams(): script to calculate position and velocity vector
% components in cylindrical coordinates at impact using Cartesian vector
% components retrieved from https://ssd.jpl.nasa.gov/sbdb.cgi 
% Arrival date: 01-26-2026 00:00:00 UTC (data available at 00:00:00)
clc,clear
day2sec = 3600*24;

% Target position
x_t = 1.386716949323424E-01; y_t = 2.722313489171105E+00;
z_t = -1.313311531032168E-01;
[th_t, r_t] = cart2pol(x_t,y_t);
r_t = vpa(r_t), th_t = vpa(th_t), z_t = vpa(z_t)
% Target velocity
vx_t = -1.059918533255927E-02/day2sec; vy_t = 1.811795239916161E-03/day2sec;
vz_t = 2.018647669674729E-04/day2sec;
vr_t = (x_t*vx_t+y_t*vy_t)/r_t; 
vth_t = (x_t*vy_t-y_t*vx_t)/r_t; 
vr_t = vpa(vr_t), vth_t = vpa(vth_t), vz_t = vpa(vz_t)

% Earth position
x_e = -5.741032104754529E-01; y_e = 7.998736632146318E-01;
z_e = -4.470917178077003E-05;
[th_e, r_e] = cart2pol(x_e,y_e);
r_e = vpa(r_e), th_e = vpa(th_e), z_e = vpa(z_e)
% Earth velocity
vx_e = -1.425646832811722E-02/day2sec; vy_e = -1.009638940775061E-02/day2sec;
vz_e = 6.804102869643572E-07/day2sec;
vr_e = (x_e*vx_e+y_e*vy_e)/r_e; 
vth_e = (x_e*vy_e-y_e*vx_e)/r_e; 
vr_e = vpa(vr_e), vth_e = vpa(vth_e), vz_e = vpa(vz_e)

% Mars position
x_m = 6.788314007926418E-01; y_m = -1.230211015379122E+00;
z_m = -4.242593103106305E-02;
[th_m, r_m] = cart2pol(x_m,y_m);
r_m = vpa(r_m), th_m = vpa(th_m), z_m = vpa(z_m)
% Mars velocity
vx_m = 1.278223105431074E-02/day2sec; vy_m = 7.962381682872869E-03/day2sec;
vz_m = -1.465771927171615E-04/day2sec;
vr_m = (x_m*vx_m+y_m*vy_m)/r_m; 
vth_m = (x_m*vy_m-y_m*vx_m)/r_m; 
vr_m = vpa(vr_m), vth_m = vpa(vth_m), vz_m = vpa(vz_m)