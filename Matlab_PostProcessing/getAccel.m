function [accelX,accelY,accelZ] = getAccel(cR,tripTime,gammaCoeff,tauCoeff,sizeC)
    gamma = angles(cR(7,1:sizeC),tripTime,gammaCoeff);
    tau = angles(cR(7,1:sizeC),tripTime,tauCoeff);
    accelR = cR(10,1:sizeC).*cos(tau).*sin(gamma);
    accelTheta = cR(10,1:sizeC).*cos(tau).*cos(gamma);
    accelZ = cR(10,1:sizeC).*sin(tau);
    accelX=accelR.*cos(cR(2,1:sizeC))-accelTheta.*sin(cR(2,1:sizeC));
    accelY=accelR.*sin(cR(2,1:sizeC))+accelTheta.*cos(cR(2,1:sizeC));

    % test Fourier calculation for start
    % vpa(gamma(1)), vpa(tau(1))
end

