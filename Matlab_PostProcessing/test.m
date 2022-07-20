clc,clear
day2sec = 3600*24;
% Asteroid position
x_a = 1.386802281980569E-01; y_a = 2.722312030081536E+00;
z_a = -1.313313155605537E-01;
[th_a, r_a] = cart2pol(x_a,y_a);
r_a = vpa(r_a), th_a = vpa(th_a), z_a = vpa(z_a)
% Asteroid velocity
vx_a = -1.059918369324184E-02/day2sec; vy_a = 1.811827137508805E-03/day2sec;
vz_a = 2.018632262875054E-04/day2sec;
vr_a = (x_a*vx_a+y_a*vy_a)/r_a; 
vth_a = (x_a*vy_a-y_a*vx_a)/r_a; 
vr_a = vpa(vr_a), vth_a = vpa(vth_a), vz_a = vpa(vz_a)

%Mars time of impact
x_m = -1.485200670893854E+00;
y_m = 7.523494901494494E-01;
z_m = 5.219918296063146E-02;
vx_m = -5.798748916827620E-03/day2sec;
vy_m = -1.128960072240110E-02/day2sec;
vz_m = -9.436678120610230E-05/day2sec;
[th_m, r_m] = cart2pol(x_m,y_m);
r_m = vpa(r_m), th_m = vpa(th_m +(2*pi)), z_m = vpa(z_m)
vr_m = (x_m*vx_m+y_m*vy_m)/r_m; 
vth_m = (x_m*vy_m-y_m*vx_m)/r_m;
vr_m = vpa(vr_m), vth_m = vpa(vth_m), vz_m = vpa(vz_m)

% Earth position
x_e = -5.740917945809353E-01; y_e = 7.998817477980121E-01;
z_e = -4.470971661342622E-05;
[th_e, r_e] = cart2pol(x_e,y_e);
r_e = vpa(r_e), th_e = vpa(th_e), z_e = vpa(z_e)
% Earth velocity
vx_e = -1.425661085314258E-02/day2sec; vy_e = -1.009619083308822E-02/day2sec;
vz_e = 6.803990745354085E-07/day2sec;
vr_e = (x_e*vx_e+y_e*vy_e)/r_e; 
vth_e = (x_e*vy_e-y_e*vx_e)/r_e; 
vr_e = vpa(vr_e), vth_e = vpa(vth_e), vz_e = vpa(vz_e)