function [] = plotData(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)

    %% Data that is required
    
    au=1.49597870691E11; % conversion of m/au
    
    % Solving differential motions
    timeFinal=(6.653820100923719e+07); % orbital period
    tspan=[timeFinal 0];
    options = odeset('RelTol',1e-12);
    [tE, yE] = ode45(@orbitalMotion,tspan,y0E,options,gammaCoeff,tauCoeff,timeFinal,0);
    [tA, yA] = ode45(@orbitalMotion,tspan,y0A,options,gammaCoeff,tauCoeff,timeFinal,0);
    
    % Transform to cartesian coordinates for position and velocity of asteroid and earth
    [cX,cY,cZ]= pol2cart(cR(2,:),cR(1,:),cR(3,:));
    [eX,eY,eZ]= pol2cart(yE(:,2),yE(:,1),yE(:,3));
    [eX_launch,eY_launch,eZ_launch]= pol2cart(launchPos(2),launchPos(1),launchPos(3));
    [aX,aY,aZ]= pol2cart(yA(:,2),yA(:,1),yA(:,3));
    % Acceleration vector in cartesian coordinates
    [accelX,accelY,accelZ] = getAccel(cR,tripTime,gammaCoeff,tauCoeff,sizeC);
    
    velX = cR(4,:).*cos(cR(2,:))-cR(1,:).*cR(5,:).*sin(cR(2,:)) - y0A(4); % AU
    velY = cR(4,:).*sin(cR(2,:))+cR(1,:).*cR(5,:).*cos(cR(2,:)) - y0A(5);
    velZ = cR(6,:) - y0A(6);
    
    %% Sub Plot 1
    
    figure(1) %orbitals
    subplot(2,3,1)
    polarplot(yE(:,2),yE(:,1),'.')
    hold on
    polarplot(yA(:,2),yA(:,1),'.')
    hold on
    polarplot(cR(2,1),cR(1,1),'r*')
    hold on
    polarplot(y0A(2),y0A(1),'*b')
    hold on
    polarplot(cR(2,:),cR(1,:))
    legend({'earth','asteroid','launch','impact','spacecraft'})
    title('r-\theta plane')
    hold off
    
    %radius vs. time
    subplot(2,3,4)
    plot(tE-(timeFinal-tripTime),yE(:,1),'.')
    hold on
    plot(tA-(timeFinal-tripTime),yA(:,1),'.')
    hold on
    plot(cR(7,:),cR(1,:))
    ylabel('r (a.u.)')
    xlabel('t (s)')
    xlim([0 tripTime])
    title('Orbital radius')
    hold off
    
    %specific angular momentum vs. time
    subplot(2,3,2)
    plot(tE-(timeFinal-tripTime),yE(:,1).*yE(:,5),'.')
    hold on
    plot(tA-(timeFinal-tripTime),yA(:,1).*yA(:,5),'.')
    hold on
    plot(cR(7,:),cR(1,:).*cR(5,:))
    ylabel('h (a.u.^{2}/s)')
    xlabel('t (s)')
    xlim([0 tripTime])
    title('Specific angular momentum')
    hold off
    
    %plot z
    subplot(2,3,3)
    plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'.')
    hold on
    plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'.')
    hold on
    plot(cR(1,:).*cos(cR(2,:)), cR(3,:),'LineWidth', 2)
    xlabel('x (a.u.)')
    ylabel('z (a.u.)')
    title('x-z plane')
    hold off
    
    %Z vs. time
    subplot(2,3,6)
    plot(tE-(timeFinal-tripTime),yE(:,3),'.')
    hold on
    plot(tA-(timeFinal-tripTime),yA(:,3),'.')
    hold on
    plot(cR(7,:),cR(3,:))
    ylabel('z (a.u.)')
    xlim([0 tripTime])
    xlabel('t (s)')
    title('Orbital elevation')
    hold off
    
    %theta vs. time
    subplot(2,3,5)
    plot(tE-(timeFinal-tripTime),mod(yE(:,2),2*pi),'.')
    hold on
    plot(tA-(timeFinal-tripTime),mod(yA(:,2), 2*pi),'.')
    hold on
    plot(cR(7,:),mod(cR(2,:), 2*pi))
    ylabel('\theta (rad.)')
    xlim([0 tripTime])
    xlabel('t (s)')
    title('Orbital angle')
    hold off
    
    %% Subplot 2
    
    figure(2) % acceleration vs. time
    subplot(2,2,1)
    plot(cR(7,:),au*cR(10,:))
    xlim([0 tripTime])
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2})')
    xlabel('t (s)')
    
    hold off
    
    %% coasting plots
    
    % figure
    % %plot(cR(7,1:sizeC),cR(11,1:sizeC))
    % xlabel('Time')
    % ylabel('coast')
    
    % figure
    % %5plot(cR(7,1:sizeC),cR(12,1:sizeC))
    % xlabel('Time')
    % ylabel('coast value')
    
    co = angles(cR(7,1:sizeC),tripTime,coast);
    coast = sin(co).^2 < coast_threshold;
    % separating the data at or above the threshold from those below
    above = sin(co).^2; below = sin(co).^2;
    above(coast) = NaN; below(~coast) = NaN;
    subplot(2,2,2)
    plot(cR(7,:),above,'color','b')
    hold on
    plot(cR(7,:),below,'color','k')
    hold on
    coast_thresholdPlot = coast_threshold*ones(1,sizeC); % creates a vector with values of coast_threshold so MATLAB can plot it as a line
    plot(cR(7,:),coast_thresholdPlot,'--','color','r')
    xlim([0 tripTime]), ylim([0,1])
    legend('thrusting','coasting','threshold')
    title('Coasting function and threshold')
    xlabel('t (s)')
    ylabel('sin^2(\psi)')
    hold off
    
    fuelSpent = (fuelMass - cR(11,:))/fuelMass;
    subplot(2,2,3)
    plot(cR(7,:),fuelSpent*100)
    xlim([0 tripTime])
    title('Fuel consumption')
    ylabel('% fuel')
    xlabel('t (s)')

    err = (cR(12,:)-cR(13,:))./cR(14,:);
    subplot(2,2,4)
    plot(cR(7,:),err*100)
    xlim([0 tripTime])
    title('Conservation of mechanical energy')
    ylabel('% error')
    xlabel('t (s)')

    % Thrust fractions and velocity components
    figure(3)

    subplot(2,3,1)
    plot(cR(7,:), au*cR(4,:))
    xlim([0 tripTime])
    title('Radial velocity')
    xlabel('t (s)')
    ylabel('v_{r} (m/s)')

    subplot(2,3,2)
    plot(cR(7,:), au*cR(5,:))
    xlim([0 tripTime])
    title('Tangential velocity')
    xlabel('t (s)')
    ylabel('v_{\theta} (m/s)')

    subplot(2,3,3)
    plot(cR(7,:), au*cR(6,:))
    xlim([0 tripTime])
    title('Axial velocity')
    xlabel('t (s)')
    ylabel('v_{z} (m/s)')
    
    subplot(2,3,4)
    plot(cR(7,:), au*cR(10,:).*sin(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime])
    %legend('matlab','c')
    title('Radial acceleration')
    xlabel('t (s)')
    ylabel('a_{r} (m/s^{2})')
    
    subplot(2,3,5)
    plot(cR(7,:), au*cR(10,:).*cos(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime])
    %legend('matlab','c')
    title('Tangential acceleration')
    xlabel('t (s)')
    ylabel('a_{\theta} (m/s^{2})')
    
    subplot(2,3,6)
    plot(cR(7,:), au*cR(10,:).*sin(cR(9,:)))
    xlim([0 tripTime])
    %legend('matlab','c')
    title('Axial acceleration')
    xlabel('t (s)')
    ylabel('a_{z} (m/s^{2})')
    
    % Thrust angle plots
    figure(4)
    
    % Test Fourier calculation for start
    % vpa(cR(8,1)), vpa(cR(9,1))
    
    subplot(2,3,1)
    plot(cR(7,:),mod(cR(8,:),2*pi))
    xlabel('t (s)'), ylabel('\gamma (rad.)')
    xlim([0 tripTime])
    title('In-plane thrust angle')
    
    subplot(2,3,2)
    plot(cR(7,:),cR(9,:))
    xlabel('t (s)'), ylabel('\tau (rad.)')
    xlim([0 tripTime])
    title('Out-of-plane thrust angle')
    
    subplot(2,3,3)
    plot(cR(7,:),co)
    xlabel('t (s)'), ylabel('\psi')
    xlim([0 tripTime])
    title('Coast series')
    
    subplot(2,3,4)
    plot(cR(7,:),sin(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime]), ylim([-1,1])
    %legend('matlab','c')
    title('Radial thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\gamma)cos(\tau)')
    
    subplot(2,3,5)
    plot(cR(7,:),cos(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime]), ylim([-1,1])
    %legend('matlab','c')
    title('Tangential thrust fraction')
    xlabel('t (s)')
    ylabel('cos(\gamma)cos(\tau)')
    
    subplot(2,3,6)
    plot(cR(7,:),sin(cR(9,:)))
    xlim([0 tripTime]), ylim([-1,1])
    %legend('matlab','c')
    title('Off-plane thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\tau)')
    
    %% full orbital plots (vectors and no vectors)
    
    radStep=1:15:length(cX)*1.0;
    %a=figure(3); % plot with vectors
    figure(5) % plot with vectors
    plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880]	)
    %plot3(cX(1),cY(1),cZ(1),'*','LineWidth', 5,'Color',[0.9290, 0.6940, 0.1250])
    %hold on
    %plot3(cX(end),cY(end),cZ(end),'*','LineWidth', 5,'Color',[0.9290, 0.6940, 0.1250])
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
    quiver3(cX(radStep),cY(radStep),cZ(radStep),accelX(radStep),accelY(radStep),accelZ(radStep),'k','Autoscalefactor',.08,'LineWidth',1)
    hold on
    [y0Ax, y0Ay, y0Az] = pol2cart(y0A(2), y0A(1), y0A(3));
    velDiff = au*sqrt((y0A(4) - cR(4,end))^2 + (y0A(5) - cR(5,end))^2 + (y0A(6) - cR(6,end))^2);
    txt = join(['t: ',num2str(tripTime/(3600*24)),' days\n|V|: ',num2str(velDiff),' m/s\nVx: ',num2str(au*velX(end)),' m/s\nVy: ',num2str(au*velY(end)),' m/s\nVz: ',num2str(au*velZ(end)),' m/s']);
    txt = compose(txt);
    text(y0Ax, y0Ay, y0Az, txt)
    title('Solar orbitals')
    legend('Spacecraft','Asteroid','Earth','Thrust')
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
    
    figure(6)
    
    r_esoi = 6.211174738e-3; % radius of Earth's sphere of influence in au
    [x,y,z] = sphere;
    
    % Earth's sphere of influence at launch
    surf(eX_launch+r_esoi*x, eY_launch+r_esoi*y, eZ_launch+r_esoi*z)
    hold on
    
    % In-plane initial position
    [alpha_x, alpha_y, alpha_z] = pol2cart(alpha, r_esoi, 0);
    plot3(alpha_x+eX_launch, alpha_y+eY_launch, alpha_z+eZ_launch,'*')
    hold on
    
    % Initial velocity vector
    quiver3(alpha_x+eX_launch, alpha_y+eY_launch, alpha_z+eZ_launch, sin(beta)*cos(zeta), cos(beta)*cos(zeta), sin(zeta),'k','Autoscalefactor',.005,'LineWidth',1);
   
    % analytical scaling
    xlim([eX_launch-2*r_esoi, eX_launch+2*r_esoi])
    ylim([eY_launch-2*r_esoi, eY_launch+2*r_esoi])
    zlim([eZ_launch-2*r_esoi, eZ_launch+2*r_esoi])
    xlabel('x (a.u.)')
    ylabel('y (a.u.)')
    zlabel('z (a.u.)')
    title('Launch conditions')
    legend({'ESOI','Position','Velocity'})
    hold off
    
    % figure (7)
    % polar(t,r_esoi*ones(1,n))
    % hold on
    % polar(alpha,r_esoi,'r*')
    % hold on
    % title('Initial position')
    % legend('ESOI','Launch')
    % hold on
    
    end