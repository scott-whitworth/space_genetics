function [] = filePlotCompare(seed1,seed2)

fileC = fopen(join(['finalOptimization-',num2str(seed1),'.bin']));
cVector = fread(fileC,Inf,'double');
fileY = fopen(join(['orbitalMotion-',num2str(seed1),'.bin']));
sizeC=cVector(end)+1;
fileD = fopen(join(['finalOptimization-',num2str(seed2),'.bin']));
dVector = fread(fileD,Inf,'double');
fileZ = fopen(join(['orbitalMotion-',num2str(seed2),'.bin']));
sizeD = dVector(end)+1;
cR = fread(fileY,[14, sizeC],'double');
dR = fread(fileZ,[14, sizeD],'double');
[tripTime1,coast_threshold1,y0E,y0T,gammaCoeff1,tauCoeff1,coast1,fuelMass1,alpha1,beta1,zeta1,launchPos1] = loadTripData(cVector);
[tripTime2,coast_threshold2,y0E,y0T,gammaCoeff2,tauCoeff2,coast2,fuelMass2,alpha2,beta2,zeta2,launchPos2] = loadTripData(dVector);
plotDataCompare(cR,y0T,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,fuelMass1,alpha1,beta1,zeta1,launchPos1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2,fuelMass2,alpha2,beta2,zeta2,launchPos2)
fclose('all');
end
