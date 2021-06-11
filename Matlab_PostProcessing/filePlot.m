function [] = filePlot(seed)

filenameO=join(['finalOptimization-',num2str(seed),'.bin']);
filenameT=join(['orbitalMotion-',num2str(seed),'.bin']);
fileC = fopen(filenameO);
cVector = fread(fileC,Inf,'double');
fileY = fopen(filenameT);
sizeC=cVector(end)+1;
cR = fread(fileY,[14, sizeC],'double');
[tripTime,coast_threshold,y0E,y0A,gammaCoeff,tauCoeff,coast,fuelMass,alpha,beta,zeta,launchPos] = loadTripData(cVector);
plotData(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
%plotDataPublish(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
fclose('all');
end