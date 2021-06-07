function yp = orbitalMotion(t,y,gammaCoeff,tauCoeff,tF,accel)

%CONSTANTS
au=1.49597870691E11; %Conversion between meters and AU
G=1.99349603314131e-44;
M=1.98892e30;

% gamma,tau, and accel as constants
% tau=3/4;
% gamma=3/4;


% gamma and tau fourier series
gamma = angles(t,tF,gammaCoeff);
tau = angles(t,tF,tauCoeff);
yp = [y(4)
    y(5)./y(1)
    y(6)
    (-G*M*y(1) ./ (y(3).^2+y(1).^2).^(3/2)) +   (y(5).^2./y(1)) + accel*cos(tau)*sin(gamma)    
    -y(4).*y(5)./y(1)                       +           0       + accel*cos(tau)*cos(gamma) 
    -G*M*y(3)./ (y(3).^2+y(1).^2).^(3/2)    +           0       + accel*sin(tau)];