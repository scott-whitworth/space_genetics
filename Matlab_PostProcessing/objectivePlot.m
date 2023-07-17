referencepoints = readtable('referencePoints-1689362410.csv');
    referencepointsx = referencepoints.Objective0;
    referencepointsy = referencepoints.Objective1;
    referencepointsz = referencepoints.Objective2;

plot3(referencepointsx,referencepointsy,referencepointsz,'o');
%plot(referencepointsx,referencepointsy,'o');

hold on

cost = readtable('1689362410-AllAdults-gen1400.csv');
    cost0 = cost.orbitPosDiffCost;
    cost1 = cost.orbitSpeedDiffCost;
    cost2 = cost.maxOrbitAssistCost;
    
    %objective0 = objective0/max(objective0);
    %objective1 = objective1/max(objective1);
    %objective2 = objective2/min(objective2);
    
plot3(cost0,cost1,cost2,'*');
%plot(cost0,cost1,'*');
xlabel('orbitPosDiffCost')
ylabel('orbitSpeedDiffCost')
zlabel('gravAssistCost')