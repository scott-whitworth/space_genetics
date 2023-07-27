function [] = plotData(cR,y0T,y0M,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
%function [] = plotData(cR,y0T,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass,alpha,beta,zeta,launchPos)
   
%% Conversion data
    %au=1.49597870691E11; % conversion of m/au
    au=1.49597870691E11; % conversion of m/au
    
    % Solving differential motions
    timeFinal=(1.5764877936e+8); % orbital period for Psyche
    %timeFinal=(3.772645011085093e+07); % orbital period for Bennu
    %timeFinal=(6.653820100923719e+07); % orbital period for Didymos
    tspan=[tripTime 0];
    %Orbial elements of Mars at t=0 (Bennu values)
    %y0M=[1.422486429200229; 0.433301411341479; -1.916394400307147E-02; 1.141483741401610E-08; 1.719277837221393E-07; 4.887178957301878E-09];
    options = odeset('RelTol',1e-13);
    [tE, yE] = ode45(@orbitalMotion,tspan,y0E,options,gammaCoeff,tauCoeff,tripTime,0);
    [tT, yT] = ode45(@orbitalMotion,tspan,y0T,options,gammaCoeff,tauCoeff,tripTime,0);
    [tM, yM] = ode45(@orbitalMotion,tspan,y0M,options,gammaCoeff,tauCoeff,tripTime,0);
    
    % Transform to cartesian coordinates for position and velocity of target and earth
    [cX,cY,cZ]= pol2cart(cR(2,:),cR(1,:),cR(3,:));
    [eX,eY,eZ]= pol2cart(yE(:,2),yE(:,1),yE(:,3));
    [eX_launch,eY_launch,eZ_launch]= pol2cart(launchPos(2),launchPos(1),launchPos(3));
    [tX,tY,tZ]= pol2cart(yT(:,2),yT(:,1),yT(:,3));
    [mX,mY,mZ]= pol2cart(yM(:,2),yM(:,1),yM(:,3));
    % Acceleration vector in cartesian coordinates
    [accelX,accelY,accelZ] = getAccel(cR,tripTime,gammaCoeff,tauCoeff,sizeC);
    
    velX = cR(4,:).*cos(cR(2,:))-cR(1,:).*cR(5,:).*sin(cR(2,:)) - y0T(4); % AU
    velY = cR(4,:).*sin(cR(2,:))+cR(1,:).*cR(5,:).*cos(cR(2,:)) - y0T(5);
    velZ = cR(6,:) - y0T(6);
      
    
    
    %% Plot 1 -------------------------------------------------------------
    
    figure(1) %orbitals
    subplot(2,3,1)
    polarplot(yE(:,2),yE(:,1),'.')
    %rlim([0 1.5])
    hold on
    polarplot(yT(:,2),yT(:,1),'.')
    hold on
    polarplot(yM(:,2),yM(:,1),'.','Color',[0.250,0.250,0.0])
    %rlim([0 1.5])
    hold on
    polarplot(cR(2,1),cR(1,1),'r*')
    %rlim([0 1.5])
    hold on
    polarplot(y0T(2),y0T(1),'*b')
    %rlim([0 1.5])
    hold on
    polarplot(cR(2,:),cR(1,:),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    %rlim([0 1.5])
    text(0, 0, '\Leftarrow Sun')
    text(cR(2,1), cR(1,1), '\leftarrow Launch')
    text(y0T(2),y0T(1), '\leftarrow Impact')
    title('r-\theta plane')
    hold off
    
    %specific angular momentum vs. time
    subplot(2,3,2)
%    plot((tE-(timeFinal-tripTime))/(3600*24),yE(:,1).*yE(:,5),'.')
    plot(tE/(3600*24),yE(:,1).*yE(:,5),'.')
    hold on
%    plot((tT-(timeFinal-tripTime))/(3600*24),yT(:,1).*yT(:,5),'.')
    plot(tT/(3600*24),yT(:,1).*yT(:,5),'.')
    hold on
    plot(tM/(3600*24),yM(:,1).*yM(:,5),'.','Color',[0.250,0.250,0.0])
    hold on
    plot(cR(7,:)/(3600*24),cR(1,:).*cR(5,:),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    ylabel('h (AU^{2}/s)')
    xlabel('t (days)')
    xlim([0 tripTime/(3600*24)])
    %ylim([1.9E-7 2.5E-7])
    legend({'Earth','target','Mars','spacecraft'},'Location','east')
    title('Specific angular momentum')
    hold off
    
    %z vs. x
    subplot(2,3,3)
    plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'.')
    hold on
    plot(yT(:,1).*cos(yT(:,2)),yT(:,3),'.')
    hold on
    plot(yM(:,1).*cos(yM(:,2)),yM(:,3),'.','Color',[0.250,0.250,0.0])
    hold on
    plot(cR(1,:).*cos(cR(2,:)), cR(3,:),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    %xlim([-3.35 3.35])
    %ylim([-0.25 0.25])
    %xlim([-2.25 2.25])
    %ylim([-0.15 0.15])
    xlabel('x (AU)')
    ylabel('z (AU)')
    title('x-z plane')
    text(cR(1,1).*cos(cR(2,1)), cR(3,1), '\leftarrow Launch')
    text(cR(1,end).*cos(cR(2,end)), cR(3,end), '\leftarrow Impact')
    hold off
    
    %radius vs. time
    subplot(2,3,4)
    plot(tE/(3600*24),yE(:,1),'.')
    hold on
    plot(tT/(3600*24),yT(:,1),'.')
    hold on
    plot(tM/(3600*24),yM(:,1),'.','Color',[0.250,0.250,0.0])
    hold on
    plot(cR(7,:)/(3600*24),cR(1,:),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    ylabel('r (AU)')
    xlabel('t (days)')
    xlim([0 tripTime/(3600*24)])
    title('Orbital radius')
    hold off
    
    %theta vs. time
    subplot(2,3,5)
    plot(tE/(3600*24),mod(yE(:,2),2*pi),'.')
    hold on
    plot(tT/(3600*24),mod(yT(:,2), 2*pi),'.')
    hold on
    plot(tM/(3600*24),mod(yM(:,2), 2*pi),'.','Color',[0.250,0.250,0.0])
    hold on
    plot(cR(7,:)/(3600*24),mod(cR(2,:), 2*pi),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    ylabel('\theta (rad)')
    xlim([0 tripTime/(3600*24)])
    xlabel('t (days)')
    title('Orbital angle')
    hold off

    %Z vs. time
    subplot(2,3,6)
    plot(tE/(3600*24),yE(:,3),'.')
    hold on
    plot(tT/(3600*24),yT(:,3),'.')
    hold on
    plot(tM/(3600*24),yM(:,3),'.','Color',[0.250,0.250,0.0])
    hold on
    plot(cR(7,:)/(3600*24),cR(3,:),'Color',[0.4660, 0.6740, 0.1880],'LineWidth', 2)
    ylabel('z (AU)')
    xlim([0 tripTime/(3600*24)])
    xlabel('t (days)')
    title('Orbital elevation')
    hold off
    
    
    %% Plot 2 -------------------------------------------------------------
    
    figure(2) % acceleration vs. time
    subplot(2,2,1)
    plot(cR(7,:)/(3600*24),au*cR(10,:),'color','k')
    xlim([0 tripTime/(3600*24)])
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2})')
    xlabel('t (days)')
    
    hold off
    
    %%coasting plots
    
    co = angles(cR(7,1:sizeC),tripTime,coast);
    coast = sin(co).^2 < coast_threshold;
    % separating the data at or above the threshold from those below
    above = sin(co).^2; below = sin(co).^2;
    above(coast) = NaN; below(~coast) = NaN;
    subplot(2,2,2)
    plot(cR(7,:)/(3600*24),above,'color','k')
    hold on
    plot(cR(7,:)/(3600*24),below,'color','k', 'LineStyle',':')
    hold on
    coast_thresholdPlot = coast_threshold*ones(1,sizeC); % creates a vector with values of coast_threshold so MATLAB can plot it as a line
    plot(cR(7,:)/(3600*24),coast_thresholdPlot,'--','color','r')
    xlim([0 tripTime/(3600*24)]), ylim([0,1])
    legend('location','southeast')
    legend('thrusting','coasting','threshold')
    title('Coasting function and threshold')
    xlabel('t (days)')
    ylabel('sin^2(\psi)')
    hold off
    
%    fuelSpent = (fuelMass - cR(11,:))/fuelMass;
%     subplot(2,2,3)
%     plot(cR(7,:)/(3600*24),fuelSpent*100,'color','k')
%     xlim([0 tripTime/(3600*24)])
%     ylim([0 100])
%     title('Fuel consumption')
%     ylabel('% fuel')
%     xlabel('t (days)')

    fuelSpent = (fuelMass - cR(11,:));
    subplot(2,2,3)
    plot(cR(7,:)/(3600*24),fuelSpent,'color','k')
    xlim([0 tripTime/(3600*24)])
    %ylim([0 100])
    title('Fuel consumption')
    ylabel('Fuel (kg)')
    xlabel('t (days)')
    hold on
    FuelSpentTxt = join(['Fuel spent=',num2str(cR(11,end),'%4.1f\n')]);
    FuelSpentTxt = compose(FuelSpentTxt);
    dim = [.2 .3 .3 .1];
    annotation('textbox',dim,'String',FuelSpentTxt,'FitBoxToText','on');
    hold off

    err = (cR(12,:)-cR(13,:))./cR(14,:);
    err(end)=0;%The last data point can have an off-by-one problem 
    subplot(2,2,4)
    plot(cR(7,:)/(3600*24),err,'color','k')
    xlim([0 tripTime/(3600*24)])
    %ylim([-0.01 0.01])
    title('Conservation of mechanical energy')
    ylabel('Error fraction')
    xlabel('t (days)')

    %% Plot 3 -------------------------------------------------------------
    % Thrust fractions and velocity components
    figure(3)

    subplot(2,3,1)
    plot(cR(7,:)/(3600*24), au*cR(4,:),'color','k')
    hold on
    plot(cR(7,end)/(3600*24),au*y0T(4,end),'*')
    xlim([0 tripTime/(3600*24)])
    title('Radial velocity')
    xlabel('t (days)')
    ylabel('v_{r} (m/s)')
    hold off

    subplot(2,3,2)
    plot(cR(7,:)/(3600*24), au*cR(5,:),'color','k')
    hold on
    plot(cR(7,end)/(3600*24),au*y0T(5,end),'*')
    xlim([0 tripTime/(3600*24)])
    title('Tangential velocity')
    xlabel('t (days)')
    ylabel('v_{\theta} (m/s)')
    hold off
    
    subplot(2,3,3)
    plot(cR(7,:)/(3600*24), au*cR(6,:),'color','k')
    hold on
    plot(cR(7,end)/(3600*24),au*y0T(6,end),'*')
    xlim([0 tripTime/(3600*24)])
    title('Axial velocity')
    xlabel('t (days)')
    ylabel('v_{z} (m/s)')
    hold off
    
    subplot(2,3,4)
    plot(cR(7,:)/(3600*24), au*cR(10,:).*sin(cR(8,:)).*cos(cR(9,:)),'color','k')
    xlim([0 tripTime/(3600*24)])
    title('Radial acceleration')
    xlabel('t (days)')
    ylabel('a_{r} (m/s^{2})')
    
    subplot(2,3,5)
    plot(cR(7,:)/(3600*24), au*cR(10,:).*cos(cR(8,:)).*cos(cR(9,:)),'color','k')
    xlim([0 tripTime/(3600*24)])
    title('Tangential acceleration')
    xlabel('t (days)')
    ylabel('a_{\theta} (m/s^{2})')
    
    subplot(2,3,6)
    plot(cR(7,:)/(3600*24), au*cR(10,:).*sin(cR(9,:)),'color','k')
    xlim([0 tripTime/(3600*24)])
    title('Axial acceleration')
    xlabel('t (days)')
    ylabel('a_{z} (m/s^{2})')
    
    %% Plot 4 -------------------------------------------------------------
    % Thrust angle plots
%     figure(4)
%     
%     % Test Fourier calculation for start
%     % vpa(cR(8,1)), vpa(cR(9,1))
%     
%     subplot(2,3,1)
%     plot(cR(7,:)/(3600*24),mod(cR(8,:),2*pi),'color','k')
%     xlabel('t (days)'), ylabel('\gamma (rad)')
%     xlim([0 tripTime/(3600*24)])
%     title('In-plane thrust angle')
%     
%     subplot(2,3,2)
%     plot(cR(7,:)/(3600*24),cR(9,:),'color','k')
%     xlabel('t (days)'), ylabel('\tau (rad)')
%     xlim([0 tripTime/(3600*24)])
%     title('Out-of-plane thrust angle')
%     
%     subplot(2,3,3)
%     plot(cR(7,:)/(3600*24),co,'color','k')
%     xlabel('t (days)'), ylabel('\psi')
%     xlim([0 tripTime/(3600*24)])
%     title('Coast series')
%     
%     subplot(2,3,4)
%     plot(cR(7,:)/(3600*24),sin(cR(8,:)).*cos(cR(9,:)),'color','k')
%     xlim([0 tripTime/(3600*24)]), ylim([-1,1])
%     title('Radial thrust fraction')
%     xlabel('t (days)')
%     ylabel('sin(\gamma)cos(\tau)')
%     
%     subplot(2,3,5)
%     plot(cR(7,:)/(3600*24),cos(cR(8,:)).*cos(cR(9,:)),'color','k')
%     xlim([0 tripTime/(3600*24)]), ylim([-1,1])
%     title('Tangential thrust fraction')
%     xlabel('t (days)')
%     ylabel('cos(\gamma)cos(\tau)')
%     
%     subplot(2,3,6)
%     plot(cR(7,:)/(3600*24),sin(cR(9,:)),'color','k')
%     xlim([0 tripTime/(3600*24)]), ylim([-1,1])
%     title('Off-plane thrust fraction')
%     xlabel('t (days)')
%     ylabel('sin(\tau)')
    
    %% Plot 5: full orbital plots (vectors and no vectors)-----------------
    
    radStep=1:15:length(cX)*1.0;
    %a=figure(3); % plot with vectors
    figure(5) % plot with vectors
    plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880])
    %xlim([-2.25 2.25])
    %ylim([-2.25 2.25])
    %zlim([-0.15 0.15])
    %xlim([-3.35 3.35])
    %ylim([-3.35 3.35])
    %zlim([-0.25 0.25])
    xlabel('x (AU)')
    ylabel('y (AU)')
    zlabel('z (AU)')
    
   
    hold on
    plot3(tX,tY,tZ,'LineWidth', 1, 'Color',	[0.6350, 0.0780, 0.1840])
    hold on
    plot3(mX,mY,mZ,'LineWidth', 1, 'Color',[0.250,0.250,0.0])
    %zlim([-2 2])
    hold on
    %plot3(eX,eY,eZ,'LineWidth', 1,'Color',[.61 .51 .74])
    plot3(eX,eY,eZ,'LineWidth', 1,'Color','b')
    %zlim([-2 2])
    hold on
    quiver3(cX(radStep),cY(radStep),cZ(radStep),accelX(radStep),accelY(radStep),accelZ(radStep),'k','Autoscalefactor',.25,'LineWidth',1)
    hold on
    [mX_final,mY_final,mZ_final]= pol2cart(y0M(2),y0M(1),y0M(3));
    r_msoi = 0.00385847; % radius of Mars' sphere of influence in AU
    [x,y,z] = sphere;
    surf(mX_final+r_msoi*x, mY_final+r_msoi*y, mZ_final+r_msoi*z, 'FaceColor','b','FaceAlpha',0.1, 'LineStyle',':')
    hold on
    [y0Tx, y0Ty, y0Tz] = pol2cart(y0T(2), y0T(1), y0T(3));

    velDiff = au*sqrt((y0T(4) - cR(4,end))^2 + (y0T(5) - cR(5,end))^2 + (y0T(6) - cR(6,end))^2);
    txt = join(['\t tripTime=',num2str(tripTime/(3600*24),'%5.2f\n'),' days\n\t |V|_{fin, rel}=',num2str(velDiff,'%5.1f'),' m/s\n\t Vx=',num2str(au*velX(end),'%5.1f'),' m/s\n\t Vy=',num2str(au*velY(end),'%5.1f'),' m/s\n\t Vz=',num2str(au*velZ(end),'%5.1f'),' m/s']);
    txt = compose(txt);
    %text(y0Tx, y0Ty, y0Tz, txt)
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',txt,'FitBoxToText','on');
    title('3-D Orbits')
    legend('Spacecraft','Target','Mars','Earth','Thrust')
    
    %r_sun = 6.96e-3; % radius of the sun in AU
    %Sun=1.0;%NOT to scale, for visualization only
    %[x,y,z] = sphere;
    %surf(Sun*r_sun*x, Sun*r_sun*y, Sun*r_sun*z, 'FaceColor','y','FaceAlpha',1.0, 'LineStyle',':')
    %legend('Sun')
    hold off

    
    %% Plot 6: Leaving ESOI -----------------------------------------------
%     figure(6)
%     
%     r_esoi = 6.211174738e-3; % radius of Earth's sphere of influence in AU
%     [x,y,z] = sphere;
%     
%     Earth=0.05;%NOT to scale, for visualization only
%     surf(eX_launch+Earth*r_esoi*x, eY_launch+Earth*r_esoi*y, eZ_launch+Earth*r_esoi*z, 'FaceColor','b','FaceAlpha',1.0, 'LineStyle',':')
%     hold on
%     % Earth's sphere of influence at launch
%     surf(eX_launch+r_esoi*x, eY_launch+r_esoi*y, eZ_launch+r_esoi*z, 'FaceColor','b','FaceAlpha',0.125, 'LineStyle',':')
%     hold on
%     
%     % In-plane initial position
%     [alpha_x, alpha_y, alpha_z] = pol2cart(alpha, r_esoi, 0);
%     plot3(alpha_x+eX_launch, alpha_y+eY_launch, alpha_z+eZ_launch,'*r')
%     hold on
%     
%     % Initial velocity vector
%     quiver3(alpha_x+eX_launch, alpha_y+eY_launch, alpha_z+eZ_launch, sin(beta)*cos(zeta), cos(beta)*cos(zeta), sin(zeta),'k','Autoscalefactor',.005,'LineWidth',1);
% 
%     % analytical scaling
%     xlim([eX_launch-2*r_esoi, eX_launch+2*r_esoi])
%     ylim([eY_launch-2*r_esoi, eY_launch+2*r_esoi])
%     zlim([eZ_launch-2*r_esoi, eZ_launch+2*r_esoi])
%     xlabel('x (AU)')
%     ylabel('y (AU)')
%     zlabel('z (AU)')
%     title('Launch conditions')
%     legend({'Earth','ESOI','Position','Velocity'})
%     
%     
%     Leave_txt = join(['alpha=',num2str(alpha,'%3.1f\n'),' rad\nbeta=',num2str(beta,'%3.1f\n'),' rad\nzeta=',num2str(zeta,'%3.1f\n'),' rad']);
%     Leave_txt = compose(Leave_txt);
%     dim = [.2 .5 .3 .3];
%     annotation('textbox',dim,'String',Leave_txt,'FitBoxToText','on');
%     
%     hold off
    
    
    end