referencepoints = readtable('referencePoints-1689707216.csv');
    referencepointsx = referencepoints.Objective0;
    referencepointsy = referencepoints.Objective1;
    referencepointsz = referencepoints.Objective2;

plot3(referencepointsx,referencepointsy,referencepointsz,'o');
%plot(referencepointsx,referencepointsy,'o');

hold on

cost = readtable('1689722765-AllAdults-gen#150.csv');
    cost0 = cost.orbitPosDiffNormalization;
    cost1 = cost.orbitSpeedDiffNormalization;
    cost2 = cost.minMarsDistNormalization;
    %cost2 = cost.maxOrbitAssistNormalization;
    
    %objective0 = objective0/max(objective0);
    %objective1 = objective1/max(objective1);
    %objective2 = objective2/min(objective2);
    
plot3(cost0,cost1,cost2,'*');
%plot(cost0,cost1,'*');
xlabel('orbitPosDiffNormalization')
ylabel('orbitSpeedDiffNormalization')
zlabel('minMarsDistNormalization')