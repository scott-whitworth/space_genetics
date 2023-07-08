function [] = filePlot(seed)

filenameO=join(['finalOptimization-',num2str(seed),'.bin']);
filenameT=join(['orbitalMotion-',num2str(seed),'.bin']);
fileC = fopen(filenameO);
cVector = fread(fileC,Inf,'double');
fileY = fopen(filenameT);
sizeC=cVector(end);
cR = fread(fileY,[14, sizeC],'double');
[tripTime,coast_threshold,y0E,y0T,y0M,gammaCoeff,tauCoeff,coast,fuelMass,alpha,beta,zeta,launchPos] = loadTripData(cVector);
plotData(cR,y0T,y0M,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
%plotDataPublish(cR,y0T,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
fclose('all');
end