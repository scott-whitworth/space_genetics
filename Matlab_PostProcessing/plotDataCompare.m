function [] = plotDataCompare(cR,y0A,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,fuelMass1,alpha1,beta1,zeta1,launchPos1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2,fuelMass2,alpha2,beta2,zeta2,launchPos2)
    %% Data that is required
    
    au=1.49597870691E11; % conversion of m/au
    
    % Solving differential motions
    timeFinal=(6.653820100923719e+07);
    tspan=[timeFinal 0];
    options = odeset('RelTol',1e-12);
    [tE, yE] = ode45(@orbitalMotion,tspan,y0E,options,gammaCoeff1,tauCoeff1,timeFinal,0);
    [tA, yA] = ode45(@orbitalMotion,tspan,y0A,options,gammaCoeff1,tauCoeff1,timeFinal,0);
    
    % Transform to cartesian coordinates for position and velocity of asteroid and earth
    [cX,cY,cZ]= pol2cart(cR(2,1:sizeC),cR(1,1:sizeC),cR(3,1:sizeC));
    [dX,dY,dZ]= pol2cart(dR(2,1:sizeD),dR(1,1:sizeD),dR(3,1:sizeD));
    [eX,eY,eZ]= pol2cart(yE(:,2),yE(:,1),yE(:,3));
    [eX_launch1,eY_launch1,eZ_launch1]= pol2cart(launchPos1(2),launchPos1(1),launchPos1(3));
    [eX_launch2,eY_launch2,eZ_launch2]= pol2cart(launchPos2(2),launchPos2(1),launchPos2(3));
    [aX,aY,aZ]= pol2cart(yA(:,2),yA(:,1),yA(:,3));
    % Acceleration vector in cartesian coordinates
    [accelX1,accelY1,accelZ1] = getAccel(cR,tripTime1,gammaCoeff1,tauCoeff1,sizeC);
    [accelX2,accelY2,accelZ2] = getAccel(dR,tripTime2,gammaCoeff2,tauCoeff2,sizeD);
    
    velX1 = cR(4,:).*cos(cR(2,:))-cR(1,:).*cR(5,:).*sin(cR(2,:)) - y0A(4); % AU
    velY1 = cR(4,:).*sin(cR(2,:))+cR(1,:).*cR(5,:).*cos(cR(2,:)) - y0A(5);
    velZ1 = cR(6,:) - y0A(6);
    velX2 = dR(4,:).*cos(dR(2,:))-dR(1,:).*dR(5,:).*sin(dR(2,:)) - y0A(4); % AU
    velY2 = dR(4,:).*sin(dR(2,:))+dR(1,:).*dR(5,:).*cos(dR(2,:)) - y0A(5);
    velZ2 = dR(6,:) - y0A(6);
    
    
    %% Sub Plot 1
    
    figure(1) %orbitals
    subplot(2,3,1)
    polarplot(yE(:,2),yE(:,1),'.')
    hold on
    polarplot(yA(:,2),yA(:,1),'.')
    hold on
    polarplot(cR(2,1),cR(1,1),'r*')
    hold on
    polarplot(dR(2,1),dR(1,1),'g*')
    hold on
    polarplot(y0A(2),y0A(1),'*b')
    hold on
    polarplot(cR(2,:),cR(1,:))
    hold on
    polarplot(dR(2,:),dR(1,:))
    legend({'earth','asteroid','launch 1','launch 2','impact','spacecraft 1','spacecraft 2'})
    title('r-\theta plane')
    hold off
    
    %radius vs. time
    minTripTime = min(tripTime1,tripTime2);
    maxTripTime = max(tripTime1,tripTime2);
    tripTimeDiff = maxTripTime - minTripTime;
    subplot(2,3,4)
    plot(tE-(timeFinal-maxTripTime),yE(:,1),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,1),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(1,:))
    else
        plot(cR(7,:),cR(1,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(1,:))
    else
        plot(dR(7,:),dR(1,:))
    end
    ylabel('r (a.u.)')
    xlabel('t (s)')
    xlim([0 maxTripTime])
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital radius')
    hold off
    
    %angular momentum vs. time
    subplot(2,3,2)
    plot(tE-(timeFinal-maxTripTime),yE(:,1).*yE(:,5),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,1).*yA(:,5),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(1,:).*cR(5,:))
    else
        plot(cR(7,:),cR(1,:).*cR(5,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(1,:).*dR(5,:))
    else
        plot(dR(7,:),dR(1,:).*dR(5,:))
    end
    ylabel('h (a.u.^{2}/s)')
    xlabel('t (s)')
    xlim([0 maxTripTime])
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Specific Angular Momentum')
    hold off
    
    %plot z
    subplot(2,3,3)
    plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'.')
    hold on
    plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'.')
    hold on
    plot(cR(1,:).*cos(cR(2,:)), cR(3,:),'LineWidth', 2)
    hold on
    plot(dR(1,:).*cos(dR(2,:)), dR(3,:),'LineWidth', 2)
    xlabel('x (a.u.)')
    ylabel('z (a.u.)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('x-z plane')
    hold off
    
    %Z vs. time
    subplot(2,3,6)
    plot(tE-(timeFinal-maxTripTime),yE(:,3),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,3),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(3,:))
    else
        plot(cR(7,:),cR(3,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(3,:))
    else
        plot(dR(7,:),dR(3,:))
    end
    ylabel('z (a.u.)')
    xlim([0 maxTripTime])
    xlabel('t (s)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital elevation')
    hold off
    
    %theta vs. time
    subplot(2,3,5)
    plot(tE-(timeFinal-maxTripTime),mod(yE(:,2),2*pi),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),mod(yA(:,2), 2*pi),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,mod(cR(2,:), 2*pi))
    else
        plot(cR(7,:),mod(cR(2,:), 2*pi))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,mod(dR(2,:), 2*pi))
    else
        plot(dR(7,:),mod(dR(2,:), 2*pi))
    end
    ylabel('\theta (rad.)')
    xlim([0 maxTripTime])
    xlabel('t (s)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital angle')
    hold off
    
    %% Subplot 2
    
    figure(2) % only spacecraft 1
    subplot(2,2,1)
    plot(cR(7,:),au*cR(10,:))
    xlim([0 tripTime1])
    legend({'spacecraft 1'})
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2})')
    xlabel('t (s)')
    
    co = angles(cR(7,:),tripTime1,coast1);
    coast1 = sin(co).^2 < coast_threshold1;
    % separating the data at or above the threshold from those below
    above1 = sin(co).^2; below1 = sin(co).^2;
    above1(coast1) = NaN; below1(~coast1) = NaN;
    subplot(2,2,2)
    plot(cR(7,:),above1,'color','b')
    hold on
    plot(cR(7,:),below1,'color','k')
    hold on
    coast_threshold1Plot = coast_threshold1*ones(1,sizeC);
    plot(cR(7,:),coast_threshold1Plot,'--','color','r')
    xlim([0 tripTime1]), ylim([0,1])
    legend({'thrusting','coasting','threshold'})
    title('Coasting function and threshold')
    xlabel('t (s)')
    ylabel('sin^2(\psi)')
    hold off
    
    fuelSpent1 = (fuelMass1 - cR(11,:))/fuelMass1;
    subplot(2,2,3)
    plot(cR(7,:),fuelSpent1*100)
    xlim([0 tripTime1])
    title('Fuel consumption')
    legend({'spacecraft 1'})
    ylabel('% fuel')
    xlabel('t (s)')

    err1 = (cR(12,:)-cR(13,:))./cR(14,:);
    subplot(2,2,4)
    plot(cR(7,:),err1*100)
    xlim([0 tripTime1])
    title('Conservation of mechanical energy')
    legend({'spacecraft 1'})
    ylabel('% error')
    xlabel('t (s)')
    
    figure(3) % only spacecraft 2
    subplot(2,2,1)
    plot(dR(7,:),au*dR(10,:))
    xlim([0 tripTime2])
    legend({'spacecraft 2'})
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2})')
    xlabel('t (s)')
    
    do = angles(dR(7,:),tripTime2,coast2);
    coast2 = sin(do).^2 < coast_threshold2;
    % separating the data at or above the threshold from those below
    above2 = sin(do).^2; below2 = sin(do).^2;
    above2(coast2) = NaN; below2(~coast2) = NaN;
    subplot(2,2,2)
    plot(dR(7,:),above2,'color','b')
    hold on
    plot(dR(7,:),below2,'color','k')
    hold on
    coast_threshold2Plot = coast_threshold2*ones(1,sizeD);
    plot(dR(7,:),coast_threshold2Plot,'--','color','r')
    xlim([0 tripTime2]), ylim([0,1])
    legend({'thrusting','coasting','threshold'})
    title('Coasting function and threshold')
    xlabel('t (s)')
    ylabel('sin^2(\psi)')
    hold off
    
    fuelSpent2 = (fuelMass2 - dR(11,:))/fuelMass2;
    subplot(2,2,3)
    plot(dR(7,:),fuelSpent2*100)
    xlim([0 tripTime2])
    title('Fuel consumption')
    legend({'spacecraft 2'})
    ylabel('% fuel')
    xlabel('t (s)')

    err2 = (dR(12,:)-dR(13,:))./dR(14,:);
    subplot(2,2,4)
    plot(dR(7,:),err2*100)
    xlim([0 tripTime2])
    title('Conservation of mechanical energy')
    legend({'spacecraft 2'})
    ylabel('% error')
    xlabel('t (s)')

    % only spacecraft 1
    figure(4)

    subplot(2,3,1)
    plot(cR(7,:), au*cR(4,:))
    xlim([0 tripTime1])
    title('Radial velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{r} (m/s)')

    subplot(2,3,2)
    plot(cR(7,:), au*cR(5,:))
    xlim([0 tripTime1])
    title('Tangential velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{\theta} (m/s)')

    subplot(2,3,3)
    plot(cR(7,:), au*cR(6,:))
    xlim([0 tripTime1])
    title('Axial velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{z} (m/s)')
    
    subplot(2,3,4)
    plot(cR(7,:), au*cR(10,:).*sin(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1])
    legend({'spacecraft 1'})
    title('Radial acceleration')
    xlabel('t (s)')
    ylabel('a_{r} (m/s^{2})')
    
    subplot(2,3,5)
    plot(cR(7,:), au*cR(10,:).*cos(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1])
    legend({'spacecraft 1'})
    title('Tangential acceleration')
    xlabel('t (s)')
    ylabel('a_{\theta} (m/s^{2})')
    
    subplot(2,3,6)
    plot(cR(7,:), au*cR(10,:).*sin(cR(9,:)))
    xlim([0 tripTime1])
    legend({'spacecraft 1'})
    title('Axial acceleration')
    xlabel('t (s)')
    ylabel('a_{z} (m/s^{2})')
    
    % only spacecraft 2
    figure(5)

    subplot(2,3,1)
    plot(dR(7,:), au*dR(4,:))
    xlim([0 tripTime2])
    title('Radial velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{r} (m/s)')

    subplot(2,3,2)
    plot(dR(7,:), au*dR(5,:))
    xlim([0 tripTime2])
    title('Tangential velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{\theta} (m/s)')

    subplot(2,3,3)
    plot(dR(7,:), au*dR(6,:))
    xlim([0 tripTime2])
    title('Axial velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{z} (m/s)')
    
    subplot(2,3,4)
    plot(dR(7,:), au*dR(10,:).*sin(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2])
    legend({'spacecraft 2'})
    title('Radial acceleration')
    xlabel('t (s)')
    ylabel('a_{r} (m/s^{2})')
    
    subplot(2,3,5)
    plot(dR(7,:), au*dR(10,:).*cos(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2])
    legend({'spacecraft 2'})
    title('Tangential acceleration')
    xlabel('t (s)')
    ylabel('a_{\theta} (m/s^{2})')
    
    subplot(2,3,6)
    plot(dR(7,:), au*dR(10,:).*sin(dR(9,:)))
    xlim([0 tripTime2])
    legend({'spacecraft 2'})
    title('Axial acceleration')
    xlabel('t (s)')
    ylabel('a_{z} (m/s^{2})')

    % Thrust angle plots (spacecraft 1)
    figure(6)
    
    % Gamma plots
    subplot(2,3,1)
    plot(cR(7,:),mod(cR(8,:),2*pi))
    xlabel('t (s)'), ylabel('\gamma (rad.)')
    xlim([0 tripTime1])
    legend('spacecraft 1')
    title('In-plane thrust angle')
    
    % Tau plots
    subplot(2,3,2)
    plot(cR(7,:),cR(9,:))
    xlabel('t (s)'), ylabel('\tau (rad.)')
    xlim([0 tripTime1])
    legend('spacecraft 1')
    title('Out-of-plane thrust angle')
    
    % Psi plots
    subplot(2,3,3)
    plot(cR(7,:),co)
    xlabel('t (s)'), ylabel('\psi')
    xlim([0 tripTime1])
    legend('spacecraft 1')
    title('Coast series')
    
    % thrust fractions

    subplot(2,3,4)
    plot(cR(7,:),sin(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Radial thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\gamma)cos(\tau)')
    
    subplot(2,3,5)
    plot(cR(7,:),cos(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Tangential thrust fraction')
    xlabel('t (s)')
    ylabel('cos(\gamma)cos(\tau)')
    
    subplot(2,3,6)
    plot(cR(7,:),sin(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Off-plane thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\tau)')

    % Thrust angle plots (spacecraft 2)
    figure(7)
    
    % Gamma plots
    subplot(2,3,1)
    plot(dR(7,:),mod(dR(8,:),2*pi))
    xlabel('t (s)'), ylabel('\gamma (rad.)')
    xlim([0 tripTime2])
    legend('spacecraft 2')
    title('In-plane thrust angle')
    
    % Tau plots
    subplot(2,3,2)
    plot(dR(7,:),dR(9,:))
    xlabel('t (s)'), ylabel('\tau (rad.)')
    xlim([0 tripTime2])
    legend('spacecraft 2')
    title('Out-of-plane thrust angle')
    
    % Psi plots
    subplot(2,3,3)
    plot(dR(7,:),do)
    xlabel('t (s)'), ylabel('\psi')
    xlim([0 tripTime2])
    legend('spacecraft 2')
    title('Coast series')
    
    % thrust fractions
    
    subplot(2,3,4)
    plot(dR(7,:),sin(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Radial thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\gamma)cos(\tau)')
    
    subplot(2,3,5)
    plot(dR(7,:),cos(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Tangential thrust fraction')
    xlabel('t (s)')
    ylabel('cos(\gamma)cos(\tau)')
    
    subplot(2,3,6)
    plot(dR(7,:),sin(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Off-plane thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\tau)')
    
    %% full orbital plots (vectors and no vectors)
    
    radStep1=1:15:length(cX)*1.0;
    radStep2=1:15:length(dX)*1.0;
    %a=figure(3); % plot with vectors
    figure(8) % plot with vectors
    plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880]	)
    hold on
    plot3(dX,dY,dZ,'LineWidth', 3,'Color',[0, 0, 1]	)
    xlim([-2.5 1.5])
    ylim([-2.0 2.0])
    zlim([-0.2 0.2])
    xlabel('x (a.u.)')
    ylabel('y (a.u.)')
    zlabel('z (a.u.)')
    
    hold on
    plot3(aX,aY,aZ,'LineWidth', 1, 'Color',	[0.6350, 0.0780, 0.1840])
    hold on
    plot3(eX,eY,eZ,'LineWidth', 1,'Color',[.61 .51 .74])
    hold on
    quiver3(cX(radStep1),cY(radStep1),cZ(radStep1),accelX1(radStep1),accelY1(radStep1),accelZ1(radStep1),'k','Autoscalefactor',.08,'LineWidth',1,'Color',[0.4660, 0.6740, 0.1880])
    hold on
    quiver3(dX(radStep2),dY(radStep2),dZ(radStep2),accelX2(radStep2),accelY2(radStep2),accelZ2(radStep2),'k','Autoscalefactor',.08,'LineWidth',1,'Color',[0, 0, 1])
    hold on
    [y0Ax, y0Ay, y0Az] = pol2cart(y0A(2), y0A(1), y0A(3));
    velDiff1 = au*sqrt((y0A(4) - cR(4,end))^2 + (y0A(5) - cR(5,end))^2 + (y0A(6) - cR(6,end))^2);
    velDiff2 = au*sqrt((y0A(4) - dR(4,end))^2 + (y0A(5) - dR(5,end))^2 + (y0A(6) - dR(6,end))^2);
    txt = join(['tripTime1: ',num2str(tripTime1/(3600*24)),' days\n|V1|: ',num2str(velDiff1),' m/s\nVx1: ',num2str(au*velX1(end)),' m/s\nVy1: ',num2str(au*velY1(end)),' m/s\nVz1: ',num2str(au*velZ1(end)),' m/s\ntripTime2: ',num2str(tripTime2/(3600*24)),' days\n|V2|: ',num2str(velDiff2),' m/s\nVx2: ',num2str(au*velX2(end)),' m/s\nVy2: ',num2str(au*velY2(end)),' m/s\nVz2: ',num2str(au*velZ2(end)),' m/s']);
    txt = compose(txt);
    text(y0Ax, y0Ay, y0Az, txt)
    title('Solar orbitals')
    legend('Spacecraft 1','Spacecraft 2','Asteroid','Earth','Thrust 1', 'Thrust 2')
    hold off
    %print(a,'3D.png','-dpng','-r300'); 
    
    
    %b=figure(4);
    %figure(4)
    %plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'LineWidth', 1,'Color',[.61 .51 .74]);
    %xlim([-2.5 2])
    %ylim([-.3 .3])
    %hold on
    %plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'LineWidth', 1,'Color',	[0.6350, 0.0780, 0.1840])
    %hold on
    %plot(cR(1,1:sizeC).*cos(cR(2,1:sizeC)), cR(3,1:sizeC),'LineWidth', 1,'Color',[0.4660, 0.6740, 0.1880])
    %xlabel('x')
    %ylabel('z')
    %legend('Earth','Asteroid','Spacecraft')
    %hold on
    %quiver(cX(radStep),cZ(radStep),accelX(radStep),accelZ(radStep),'k','LineWidth',1,'MarkerSize',15,'Autoscalefactor',.08)
    %title('Solar orbitals')
    %hold off
    %print(b,'2DNoVec.png','-dpng','-r350'); 

    figure(9)
    
    r_esoi = 6.211174738e-3; % radius of Earth's sphere of influence in au
    [x,y,z] = sphere;
    
    % spacecraft 1
    
    % Earth's sphere of influence at launch
    surf(eX_launch1+r_esoi*x, eY_launch1+r_esoi*y, eZ_launch1+r_esoi*z)
    hold on
    % In-plane initial position
    [alpha_x1, alpha_y1, alpha_z1] = pol2cart(alpha1, r_esoi, 0);
    plot3(alpha_x1+eX_launch1, alpha_y1+eY_launch1, alpha_z1+eZ_launch1,'*','Color',[0 0.5 0.5])
    hold on
    % Initial velocity vector
    quiver3(alpha_x1+eX_launch1, alpha_y1+eY_launch1, alpha_z1+eZ_launch1, sin(beta1)*cos(zeta1), cos(beta1)*cos(zeta1), sin(zeta1),'Autoscalefactor',.005,'LineWidth',1,'Color',[0 0.5 0.5]);
    hold on
    
    % spacecraft 2
    
    % Earth's sphere of influence at launch
    surf(eX_launch2+r_esoi*x, eY_launch2+r_esoi*y, eZ_launch2+r_esoi*z)
    hold on
    % In-plane initial position
    [alpha_x2, alpha_y2, alpha_z2] = pol2cart(alpha2, r_esoi, 0);
    plot3(alpha_x2+eX_launch2, alpha_y2+eY_launch2, alpha_z2+eZ_launch2,'*','Color',[0.5 0 0.5])
    hold on
    % Initial velocity vector
    quiver3(alpha_x2+eX_launch2, alpha_y2+eY_launch2, alpha_z2+eZ_launch2, sin(beta2)*cos(zeta2), cos(beta2)*cos(zeta2), sin(zeta2),'Autoscalefactor',.005,'LineWidth',1,'Color',[0.5 0 0.5]);
   
    % analytical scaling
    max_x = eX_launch1; max_y = eY_launch1; max_z = eZ_launch1;
    min_x = eX_launch2; min_y = eY_launch2; min_z = eZ_launch2;
    if max_x < min_x
        max_x = eX_launch2;
        min_x = eX_launch1;
    end
    if max_y < min_y
        max_y = eY_launch2;
        min_y = eY_launch1;
    end
    if max_z < min_z
        max_z = eZ_launch2;
        min_z = eZ_launch1;
    end
    diff_x = max_x - min_x; diff_y = max_y - min_y; diff_z = max_z - min_z;
    if diff_x > diff_y && diff_x > diff_z
        diff_max = diff_x;
    elseif diff_y > diff_z
        diff_max = diff_y;
    else
        diff_max = diff_z;
    end
    scale_x = (diff_max-diff_x)/2 + 2*r_esoi;
    scale_y = (diff_max-diff_y)/2 + 2*r_esoi; 
    scale_z = (diff_max-diff_z)/2 + 2*r_esoi;
    xlim([min_x-scale_x, max_x+scale_x])
    ylim([min_y-scale_y, max_y+scale_y])
    zlim([min_z-scale_z, max_z+scale_z])
    % % Test limits
    % xlim([min_x-2*r_esoi, max_x+2*r_esoi])
    % ylim([min_y-2*r_esoi, max_y+2*r_esoi])
    % zlim([min_z-2*r_esoi, max_z+2*r_esoi]) 
    
    xlabel('x (a.u.)')
    ylabel('y (a.u.)')
    zlabel('z (a.u.)')
    title('Launch conditions')
    legend({'ESOI 1','Position 1','Velocity 1','ESOI 2','Position 2','Velocity 2'})
    hold off
    
    end
