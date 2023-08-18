referencepoints = readtable('referencePoints-1692231162.csv');
    referencepointsx = referencepoints.Objective0;
    referencepointsy = referencepoints.Objective1;
    referencepointsz = referencepoints.Objective2;

plot3(referencepointsx,referencepointsy,referencepointsz,'o');
%plot(referencepointsx,referencepointsy,'o');

hold on

cost = readtable('1692231162-AllAdults-gen#200.csv');
    cost0 = cost.posDiffNormalization;
    cost1 = cost.speedDiffNormalization;
    cost2 = cost.fuel_spentNormalization;
    %cost2 = cost.maxOrbitAssistNormalization;
    
    %objective0 = objective0/max(objective0);
    %objective1 = objective1/max(objective1);
    %objective2 = objective2/min(objective2);
    
plot3(cost0,cost1,cost2,'*');
%plot(cost0,cost1,'*');
xlabel('posDiffNormalization')
ylabel('speedDiffNormalization')
zlabel('fuel_spentNormalization')