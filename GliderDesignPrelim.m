%% Housekeeping
clc 
clear all
close all

%% Provided flat plate CL and CD values with respect to AOA
plateCl = [-.75 -.76 -.77 -.78 -0.8 -0.78 -.76 -.74 -.7 -.6 -.5 -.38 -.3 -.2 -.1 0 .1 .2 .3 .4 .49 .6 .7 .75 .76 .77 .76 .75 .75 .74 .74 .73 .72 .71 .71 .7];
plateCd = [.21 .2 .19 .18 .17 .15 .13 .11 .09 .07 .05 .04 .02 .01 .005 0  .005 .01 .015 .02 .05 .065 0.09 .11 .13 .15 .17 .19 .199 .21 .22 .23 .25 .27 .285 .29];
plateAOA = [-15:1:20];
plateAOAMod = plateAOA(plateAOA >= 0 & plateAOA <= 12);
plateClMod = plateCl(plateAOA >= 0 & plateAOA <= 12);

%% Constants
AR = 6;

g = 9.81;
rho_air = 1.0581; %kg / m^3
rho_foam = .295*g; %N/m^2

Swetfus = 1115.784/(100^2); %meters^2     1115.784/(100^2);
Swetfus = Swetfus + .33;
Swettail = ((0.35* 0.18) + (.2*0.08)) * 2;
Swetrest = Swetfus + Swettail; 

W_cam = .160 * g; %weight in N
W_fus = .2791 + 0.33 /2 *rho_foam; %N 
W_tail = Swettail/2 * rho_foam; %N
W_B =  0.1*9.81;
W_rest = W_fus + W_tail + W_cam + W_B; %total

cfe = 0.0055;
e = 1.46;
eo = (1.78*(1-(0.045*(AR^0.68)))) - 0.64;
k = 1 / (pi * eo *AR);
h0 = 7; %initial height
meu = 1.74*10^(-5)
h = 0.06; %m

S_refMax = 4/AR; 

% all % 
S_ref = 0.01 : 0.01 : S_refMax;
b = sqrt(S_ref.*AR);
W = W_rest + S_ref*rho_foam; 
WL = W./S_ref; 
CD0 = cfe.*((2.*S_ref + Swetrest)./S_ref);

D_ratio = (0.44+(0.9594*h./(b)))./(0.44+(2.219*h./(b)));

%% Range
%range

CL_r = sqrt(CD0./(D_ratio.*k));
CD_r = 2*CD0;

R = (CL_r./CD_r)*h0; 
[a, idx] = min(abs(R-100))
R_ideal = R(idx)
v_r = getvelo(CL_r, rho_air, W, S_ref); 

%% Endurance 
CL_e = sqrt(3*CD0./(D_ratio.*k));
CD_e = 4*CD0;
v_e = getvelo(CL_e, rho_air, W, S_ref); 
D_e = getdrag(CD_e, rho_air, v_e, S_ref);
SinkRate = v_e.*(D_e./W);
Endurance = h0 ./ SinkRate; 


%% Estimating Lift Curve 
a_0 = 0.1; %2-D CL slope
CLalpha = a_0 / (1+((57.3*a_0)/(pi*e*AR))); 
alpha = [0:1:9]; 
CL = 0.5*alpha.*CLalpha + 0.5*alpha.*CLalpha*0.9; %stalls @ AoA = 9
CL_max = max(CL);
stall_v = getvelo(CL_max,rho_air,W,S_ref);

%% Drag Polar and Final Design
S_ref_final = S_ref(idx)
b_final = b(idx)
c_final = S_ref_final/b_final

CD = CD0(idx) + D_ratio(idx) * k * CL.^2;

alpha_final = CL_r(idx)/CLalpha
%% Plotting
%Rmax, Rmin, Rmax VS WL, CL vs WL
figure(1) 
plot(WL, R, 'b')
hold on
yline(100, 'r--')
yline(70, 'r--')
ylabel('Range (m)')
xlabel('Wing Loading Required (N/m^2)')
yyaxis right
ylabel ('C_L Required')
plot(WL, CL_r, 'g')
yline(CL_max)
legend('Range', 'Max Range Req', 'Min Range Req','C_L Required', 'CL max')
title('Max R & CL req vs W/S')
xlim([0 45])
saveas(gcf, 'Graphs/MaxR_CLReq.png')

%Emax req, Emin req, E CL req for E, CL max limit
figure(2)
plot(WL, Endurance, 'b')
hold on 
yline(10, '--r')
yline(7, '--r') 
ylabel('Endurance Acheived (s)')
yyaxis right
ylabel ('C_L Required')
plot(WL, CL_e, 'g') 
yline(CL_max)
xlabel('Wing Loading Required (N/m^2)')
title('Max E & CL req vs W/S')
legend('Endurance', 'Max Endurance Req', 'Min Endurance Req','C_L Required', 'CL max')
xlim([0, 45])
saveas(gcf, 'Graphs/MaxE_CLReq.png')

%Vel Req for Max R & Max E
figure(3)
plot(WL, v_r, 'b') %velocity for range
hold on 
plot(WL, v_e, 'g') %velocity for endurance 
yline(12, '--r')
yline(7, '-.r')
plot(WL, stall_v, '--b')
ylabel('Velocity Required (m/s)')
xlabel('Wing Loading Required (N/m^2)')
title('Vel Req for Max R & Max E')
legend('Velo for Rmax', 'Velo for Emax', 'Max Velo Req','Min Velo Req', 'Stall velo for CL max')
xlim([0, 45])
saveas(gcf, 'Graphs/VelMaxEMaxRReq.png')

%Variation of aircraft weight with wing loading 
figure(4)
plot(WL, W)
xlabel('Wing Loading Required (N/m^2)')
ylabel('Aircraft Weight (N)') 
title('Variation of Aircraft Weight with W/S')
xlim([0, 45])
saveas(gcf, 'Graphs/AircraftWeight_WS.png')

%CL vs CD
figure(5)
plot(CL,CD)
xlabel('CL')
ylabel('CD')
title('CD vs CL')
saveas(gcf, 'Graphs/CDvsCL.png')

%AoA vs CL
figure(6)
plot(plateAOAMod, plateClMod)
hold on;
plot(alpha,CL)
hold off;
xlabel('AoA (degrees)')
ylabel('CL') 
title('CL vs Angle of Attack')
legend('2D Flat Plate', '3D Wing Design', 'location', 'northwest')
saveas(gcf, 'Graphs/CLvsAOA.png')

%Variation of aircraft cost with wing loading 
figure(7)
plot(WL, (W/g)*1000)
xlabel('Wing Loading Required (N/m^2)')
ylabel('Aircraft Cost ($)') 
title('Variation of Aircraft Cost with W/S')
xlim([0, 45])
saveas(gcf, 'Graphs/AircraftCost_WS.png')

%Geometry
figure(8)
plot(WL, b/2,'linewidth', 1.25) %wingspan vs WL
hold on
yline(1) %Wingspan limit
plot(WL, getChord(S_ref, b), '--','linewidth', 1.25) % chord length, no taper - root = tip
plot(WL, getChord(S_ref, b)./4, '-.','linewidth', 1.25) %aerodynamic center distance from leading edge (1/4 chord) 
plot(WL, getXcr(v_r, rho_air, meu), '--','linewidth', 1.25) %Xcr
xlabel('Wing Loading Required (N/m^2)')
ylabel('Length (m)') 
ylim([0, 1.5])
yyaxis right
plot(WL, S_ref,'linewidth', 1.25) %Planform area vs WL
ylabel('Planform Area - Sref (m^2)')
title('Variation of Wing Geometry with Constant Taper and AR')
legend('Wing Span', 'Wing Span Limit', 'Chord Length (No Taper)', 'MAC', 'Laminar Trans Xcr', 'Planform Area') 
xlim([0, 45])
saveas(gcf, 'Graphs/WingGeoConstTaperAR.png')

%% Function 
%calculates laminar transition point
function out = getXcr(v, rho, meu)
Re_cr = 5*10^5; 
out = ( Re_cr .* meu ) ./ ( rho .* v );
end

%gets chord length
function out = getChord(S_ref, b)
out = S_ref./b;
end

% gets velocity using L = W
function out = getvelo(CL,rho,W,S)
out = sqrt(W./((1/2).*rho.*S.*CL));
end

% gets velocity using L = W
function out = getdrag(CD,rho,v,S)
out = 1/2.*CD.*v.^2.*rho.*S;
end