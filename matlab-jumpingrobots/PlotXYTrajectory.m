function [h,d,velocity] = PlotXYTrajectory(energy, mass, rho, Cd, A, angle, g)
% PlotXYTrajectory(energy, mass, rho, Cd, A, angle, g)
%
% Primarily to test ODE function to plot a trajectory for the given values.
% energy = vector of energies in J
% mass = mass in kg
% rho = density of air in kg/m^3
% Cd = drag coefficient
% A = cross-sectional area in m^2
% angle = takeoff angle in degrees
% g = gravity in m/s^2

c='bgrkmyc';
s(1,:) = '- ';
s(2,:) = ': ';
s(3,:) = '--';
s(4,:) = '-.';
s(5,:) = '+ ';
s(6,:) = 'o ';
s(7,:) = 's ';
clf;

legendstrings = cellstr(int2str(zeros(2*length(energy),1)));
for i=1:length(energy)

    v = sqrt(2*energy(i)/mass);
    viy = v*sin(deg2rad(angle));
    vix = v*cos(deg2rad(angle));
    
    tf = 2*viy/g;
    [t,x] = ode45(@TrajectoryWithoutDrag,[0,tf],[0,0,vix,viy],[],g);
    [tdrag,xdrag] = ode45(@TrajectoryWithXYDrag,[0,tf],[0,0,vix,viy],[],mass,rho*Cd*A/2,g);
%     [tdrag,xdrag] = ode45(@TrajectoryWithXYAndViscousDrag,[0,tf],[0,0,vix,viy],[],mass,rho*Cd*A/2,sqrt(A),g);
        
    x1 = x(:,1); y1 = x(:,2);
    x2 = xdrag(:,1);
    y2 = xdrag(:,2);
    
    beta = Cd*rho*A/2;

    tfdrag = 2*sqrt(mass/beta/g)*atan(sqrt(beta/g/mass)*viy);
    tdrag2 = 0:1e-5:tfdrag;
    
    a = sqrt(beta*g/mass);
    b = sqrt(beta/g/mass);
    x3 = (mass/beta).*log(1+(beta/mass)*vix.*tdrag2);
    y3 = (mass/beta).*log(cos(a*tdrag2)+b*viy*sin(a*tdrag2));

    vxdrag = xdrag(:,3); vydrag = xdrag(:,4);
    % Modify xdrag due to tf being longer than needed here - remove any
    % pieces from where xdrag(:,2) < 0
    indicesLessThanZero = find(y2 < 0);
    if size(indicesLessThanZero) ~= 0
        tdrag = tdrag(1:(min(indicesLessThanZero)));
        y2 = y2(1:(min(indicesLessThanZero)));
        x2 = x2(1:(min(indicesLessThanZero)));
        vxdragf = vxdrag(min(indicesLessThanZero));
        vydragf = vydrag(min(indicesLessThanZero));
    else
        vxdragf = vxdrag(end);
        vydragf = vydrag(end);        
    end
    vdragf = sqrt(vxdragf^2 + vydragf^2);
    
    [hvac, ihvac] = max(y1);
    dvac = x1(ihvac);
    [hair, ihair] = max(y2);
    dair = x2(ihair);
    h(i,:) = [hvac*100,hair*100];
    d(i,:) = [max(x1)*100,max(x2)*100];
    velocity(i,:) = [v,vdragf];
    d3 = 100*(mass/beta)*log(1+2*vix*sqrt(beta/mass/g)*atan(sqrt(beta/mass/g)*viy));
    h3 = 100*(mass/beta/2)*log(1+(beta/mass/g)*viy^2);
    
    if mod(i,7) == 0
        color = c(7);
    else
        color = c(mod(i,7));
    end
    if mod(i,7) == 0
        style = s(7,:);
    else
        style = s(mod(i,7),:);
    end
    hold on;
    plothandles(2*i-1) = plot(x1*100,y1*100,strcat(color,style),'linewidth',5);
    plothandles(2*i) = plot(x2*100,y2*100,strcat(color),'linewidth',5);
%     plot(x3*100,y3*100,strcat(color),'linewidth',5);
    legendstrings(2*i-1) = cellstr(strcat('U_{kinetic} = ',int2str(energy(i)*1e6),' \muJ in vacuum'));
    legendstrings(2*i) = cellstr(strcat('U_{kinetic} = ',int2str(energy(i)*1e6),' \muJ in air'));
    plot(dair*100,hair*100,strcat(color,'x'),'linewidth',20)
%     vactext = text(dvac*100,hvac*100+.25,sprintf('h = %0.1f cm\nd = %0.1f cm',hvac*100,max(x1)*100));
%     set(vactext,'FontSize',18)
    airtext = text(dair*100+1,hair*100+.6,sprintf('h = %0.1f cm\nd = %0.1f cm',h3,d3));
    set(airtext,'FontSize',18)
    
end

newaxis = axis;
newaxis(3) = 0;
axis(newaxis)

titlestring = strcat(sprintf('Jumping Trajectory, Mass = %0.3g mg, Angle = %0.3g',mass*1e6,angle),'^\circ');
title(titlestring,'FontSize',24);
xlabel('distance (cm)','FontSize',24);
ylabel('height (cm)','FontSize',24);

legend(plothandles(:),legendstrings,0);
set(gca,'FontSize',18)

hold off;
