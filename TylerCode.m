clear
clc

% %% Equations
% rho = 1.0581; % density at alt 1.5 km (kg/m^3)
% AR = nan;
% rhoM = 0.295; % g/m^3
% g = 9.81;
% e = 0.9;
% e0 = (1.78*(1-(0.045*(AR^0.68)))) - 0.64;
% h = 7;
% 
% % Weights
% syms Sref
% Wfuse = nan;
% WHTail = nan;
% WVTail = nan;
% WBallast = nan;
% WPay = nan;
% Wwing = Sref*rhoM*g;
% % 2
% Wrest = Wfuse + WHTail + WVTail + WBallast + WPay;
% WTO = Wfuse + Wwing + WHTail + WVTail + WBallast + WPay;
% % WTO = Wrest + Sref*rhoM*g
% % WingLoad = WTO/Sref = Wrest/Sref+rhoM*g
% 
% % Lift
% alpha = -10:20;
% a0 = nan;
% alpha0 = 0;
% a = a0/(1+a0/(pi*AR*e));
% CL = a*(alpha-alpha0);
% 
% % Drag
% Swet = nan;
% Cfe = 0.0030;
% CDmin = Cfe*(Swet/Sref);
% % 1
% CD0 = CDmin;
% CD = CD0 + CL.^2/(pi*e0*AR);
% % CD = Cfe*(Swet/Sref) + CL(alpha)^2/(pi*e0*AR)
% 
% 
% 
% %% Performance Graphs
% syms WingLoad
% % Range
% % Rmax = h*LDmax
% % Rmax = h*CL(AOALDmaxIndex)/(Cfe*(Swet/Sref) + CL(AOALDmaxIndex)^2/(pi*e0*AR))
% Rmax = h*CL(AOALDmaxIndex)/(WingLoad*Cfe*Swet/WTO + CL(AOALDmaxIndex)^2/(pi*e0*AR));
% % vRmax = sqrt(2*WTO/(CL(AOALDmaxIndex)*Sref*rho))
% % vRmax = sqrt((WTO/Sref)*2/(CL(AOALDmaxIndex)*rho))
% vRmax = sqrt(WingLoad*2/(CL(AOALDmaxIndex)*rho));
% % Endurance
% % vEmax = sqrt(2*WTO/(CL(AOAPRminIndex)*Sref*rho))
% % vEmax = sqrt(WTO/Sref*2/(CL(AOAPRminIndex)*rho))
% vEmax = sqrt(WingLoad*2/(CL(AOAPRminIndex)*rho));
% % DEmax = CDPRmin*1/2*vEmax^2*rho*Sref
% % DEmax = (Cfe*(Swet/Sref) + CL(AOAPRminIndex)^2/(pi*e0*AR))*1/2*(WTO/Sref*2/(CL(AOAPRminIndex)*rho))*rho*Sref
% % DEmax = (Cfe*(Swet/Sref) + CL(AOAPRminIndex)^2/(pi*e0*AR))*(WTO/(CL(AOAPRminIndex)))
% DEmax = (WingLoad*Cfe*Swet/CL(AOAPRminIndex) + WTO*CL(AOAPRminIndex))/(pi*e0*AR);
% % sink = vEmax*DEmax/WTO
% % sink = sqrt(WingLoad*2/(CL(AOALDmaxIndex)*rho))*((WingLoad*Cfe*Swet/CL(AOAPRminIndex) + WTO*CL(AOAPRminIndex))/(pi*e0*AR))/WTO
% sink = sqrt(WingLoad*2/(CL(AOALDmaxIndex)*rho))*((1/Sref*Cfe*Swet/CL(AOAPRminIndex) + CL(AOAPRminIndex))/(pi*e0*AR));
% Emax = h/sink;

% Initializations
h = 7;
% b = 1; % (max) wingspan
% AR = b^2/Sref; % contsant!
AR = 10;
e0 = 1;
Cfe = 0.0030;
k = 1/(pi*AR*e0);
rhoM = 0.295; % density of foam (N/m^2)???
rho = 1.0581; % density at alt 1.5 km (kg/m^3)
Wrest = 0.5; % weight of glider minus wings (N)
% Swet = 2*Sref+Srest;
Srest = 0.5; % surafce area of fuselage and tail

% Range: R = (L/D)*h, L/Dmax at CD0 = kCl^2
syms Sref
SrefVec = 0.01:0.01:1;
CD0 = Cfe*(2+Srest/Sref);
CL_LDmax = sqrt(CD0/k);
CL_LDmaxPlot = double(subs(CL_LDmax,Sref,SrefVec));
CD_LDmax = 2*CD0;
WTOI = rhoM*Sref + Wrest; % syms WTO
v_LDmax = sqrt(2*WTOI/(CL_LDmax*rho*Sref));
vRmaxPlot = double(subs(v_LDmax,Sref,SrefVec));
Rmax = (CL_LDmax/CD_LDmax)*h;
RmaxPlot = double(subs(Rmax,Sref,SrefVec)); % max range in terms of Sref
% WTO = rhoM*SrefVec + Wrest; % total weight
WTO = double(subs(WTOI,Sref,SrefVec)); % total weight
RmaxWL = WTO./SrefVec; % wing loading for given Sref

% Endurance: E = h/sink, sink = v*D/W, PRmin: CD0 = 1/3*k*CL^2
syms Sref
CD0 = Cfe*(2+Srest/Sref);
CL_PRmin = sqrt(3*CD0/k);
CL_PRminPlot = double(subs(CL_PRmin,Sref,SrefVec));
CD_PRmin = CD0*4/3;
WTOI = rhoM*Sref + Wrest; % syms WTO
v_PRmin = sqrt(2*WTOI/(CL_PRmin*rho*Sref));
vEmaxPlot = double(subs(v_PRmin,Sref,SrefVec));
D_PRmin = CD_PRmin*1/2*rho*v_PRmin.^2*Sref;
sink = v_PRmin*D_PRmin/WTOI;
Emax = h/sink;
EmaxPlot = double(subs(Emax,Sref,SrefVec));
EmaxWL = WTO./SrefVec;

% Finding stall AoA and CL limit
e = 0.9;
a0 = 0.1;
alpha0 = 0;
a = a0/(1+a0/(pi*AR*e));
% CL = a*(alpha-alpha0);
alphaStall = 8; % (rad), 8 or 9 degrees? 
CLlim = a*(alphaStall-alpha0);

figure(1)
plot(RmaxWL,RmaxPlot)
hold on
plot(RmaxWL,70*ones(1,length(RmaxWL)))
plot(RmaxWL,100*ones(1,length(RmaxWL)))
hold off
ylabel('Range (m)')
yyaxis right
plot(RmaxWL,CL_LDmaxPlot)
hold on
plot(RmaxWL,CLlim*ones(1,length(RmaxWL)))
hold off
xlabel('Wing Loading (N/m^2)')
ylabel('Coefficient of Lift Required')
title('Maximum Range vs Wing Loading')
legend('Range Achieved','Min Range Req','Max Range Req','CL Req for Range','CL Max Limit')

figure(2)
plot(EmaxWL,EmaxPlot)
hold on
plot(EmaxWL,7*ones(1,length(EmaxWL)))
plot(EmaxWL,10*ones(1,length(EmaxWL)))
ylabel('Endurance (s)')
yyaxis right
plot(EmaxWL,CL_PRminPlot)
plot(EmaxWL,CLlim*ones(1,length(EmaxWL)))
hold off
title('Maximum Endurance vs Wing Loading')
xlabel('Wing Loading (N/m^2)')
ylabel('Coefficient of Lift Required')
legend('Endurance Achieved','Min Endurance Req',' Max Endurance Req','CL Req for Endurance','CL Max Limit')

figure(3)
plot(EmaxWL,vEmaxPlot)
hold on
plot(RmaxWL,vRmaxPlot)
hold off
title('Velocity Required for Maximum Range and Endurance')
xlabel('Wing Loading (N/m^2)')
ylabel('Required Velocity (m/s)')
legend('Velocity for Maximum Endurance','Velocity for Maximum Range')

figure(4)
plot(RmaxWL,WTO)
title('Glider Weight vs Wing Loading')
xlabel('Wing Loading (N/m^2)')
ylabel('Glider Weight (N)')

figure(5)
% plot(RmaxWL,cost)
title('Glider Cost vs Wing Loading')
xlabel('Wing Loading (N/m^2)')
ylabel('Glider Cost ($)')

figure(6)
title('Glider Wing Geometry vs Wing Loading')
xlabel('Wing Loading (N/m^2)')
ylabel('Length (m)')
yyaxis right
plot(RmaxWL,SrefVec)
ylabel('Planform Area (m^2)')