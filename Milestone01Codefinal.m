%% Lab 01 - Milestone 1 - Tempest UA

%% Housekeeping
clear
clc 
close all
load("TempestWorkSpace") 
%% Data
%constants 
Re = 200000;            %reynolds #
e = 0.90;               %assumed efficiency factor
AR = 16.5;              %aspect ratio
mu = 1.74*(10^(-5));    %viscosity at alt 1.5 km (kg/m s)
rho = 1.0581;           %density at alt 1.5 km (kg/m^3)
v = 21;
GTOW = 6.4; %kg
S = 0.63; %m^2

%AoAs for airfoildata and drag polar data
alpha_AF = MH32AFdata{1:18,1};
alpha_DP = MH32AFdata{1:18,1};
%cl and cd for airfoil data
Cl = MH32AFdata{1:18,2};
Cd = MH32AFdata{1:18,3};
%CL and CD for drag polar data 2-D
CL_true = DragPolarData{1:18,2};
CD_true = DragPolarData{1:18,3};

%% 2-D to 3-D approximation (part 1)

Cl_fit = polyfit(alpha_AF(6:14), Cl(6:14), 1);
a0 = Cl_fit(1);                             %lift curve slope (C_L vs alpha) 2D
a = a0 / (1 + ((57.3*a0) / (pi * e * AR))); %lift curve slope 3D approx

[~, idx] = min(abs(Cl));                    %finds index where lift is nearest zero 
alpha_zero = alpha_AF(idx);                 %AoA where lift is zero

CL = a * (alpha_AF - alpha_zero);           %CL 3-D approximation 

%% Drag Polars (part 2)
CDwing = Cd + (CL.^2 / (pi * e * AR));    %3D wing drag polar - calculated by adding induced drag to 2D profile drag (Cd)

e0 = (1.78*(1-(0.045*(AR^0.68)))) - 0.64;  
k1 = 1/(pi*e0*AR);
Cfe = 0.0055;
Swet = 2.7306;
Sref = 0.63; %clarify with TA if correct

CDmin = Cfe*(Swet/Sref);            %(equation 10)
[~, i] = min(CDwing);
CLminD = CL(i);                     %CL corresponding to CDmin

CD_0 = CDmin + (k1*(CLminD.^2));    %parasite drag (equation 9)
CDwhole6 = CDmin + (k1*((CL-CLminD).^2));               %whole aircraft drag polar (equation 6) 

tail = 2*altF(0.13, 0.75, v, mu, rho);
wing = 4*altF(0.23, 1.53, v, mu, rho);
vert = 2*altF(0.2, 0.4, v, mu, rho);
body = 2*altF(1.56, 2*pi*0.08, v, mu, rho);

CDminalt = tail+wing+vert+body;
CDwhole6alt = CDminalt + (k1*((CL-CLminD).^2));



%% Performance
%% approx
%part i. 
LDMcalc = CL./CDwhole6; %L/D 
[LDMcalc, LDMcalci] = max(LDMcalc)
% LDMcalc_fit = polyfit()
glideR_v_calc = getvelo(CL(LDMcalci),rho,GTOW,S)

%part ii. 
%Max prop powered range also occurs at L/D max
powerR_v_calc = getvelo(CL(LDMcalci),rho,GTOW,S)

%part iii
LDMEcalc = (CL).^(3/2)./CDwhole6; %L/D 
[LDMEcalc, LDMEcalci] = max(LDMEcalc);
powerE_v_calc = getvelo(CL(LDMEcalci),rho,GTOW,S)

%% flat plate approx
%part i. 
LDMcalcalt = CL./CDwhole6alt; %L/D 
[LDMcalcalt, LDMcalcialt] = max(LDMcalcalt)
% LDMcalc_fit = polyfit()
glideR_v_calcalt = getvelo(CL(LDMcalcialt),rho,GTOW,S)

%part ii. 
%Max prop powered range also occurs at L/D max
powerR_v_calcalt = getvelo(CL(LDMcalcialt),rho,GTOW,S)

%part iii
LDMEcalcalt = (CL).^(3/2)./CDwhole6alt; %L/D 
[LDMEcalcalt, LDMEcalcialt] = max(LDMEcalcalt);
powerE_v_calcalt = getvelo(CL(LDMEcalcialt),rho,GTOW,S)




%% true 
%part i.
LDMtrue = CL_true./CD_true; %L/D 
[LDMtrue, LDMtruei] = max(LDMtrue)
glideR_v_true = getvelo(CL_true(LDMtruei),rho,GTOW,S)

%part ii. 
%Max prop powered range also occurs at L/D max
powerR_v_true = getvelo(CL_true(LDMtruei),rho,GTOW,S)

%part iii
LDMEtrue = (CL_true).^(3/2)./CD_true; %L/D 
[LDMEtrue, LDMEtruei] = max(LDMEtrue);
powerE_v_true = getvelo(CL(LDMEtruei),rho,GTOW,S)



%% Plotting

figure(1) %CL vs alpha for 2D and 3D 
plot(alpha_AF, Cl, "r")   %2D 
grid on
yline(0);
xline(0); 
hold on
plot(alpha_AF, CL, "b")   %3D
hold off
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift')
legend('2D airfoil','3D approximation')
title('Tempest Coefficient of Lift vs. Angle of Attack')


 figure(2) %CD vs CL 
plot(CL_true, CD_true, "g")
hold on
plot(CL, CDwing, "k")
plot(CL, CDwhole6, "b")
% xlabel("CL")
% ylabel("CD")
% yline(0)
% xline(0)
% legend("True", "Finite Wing", "Whole AC Approx")
% title('Tempest CL vs. CD')

%figure(3) %CD vs CL alt
% plot(CL_true, CD_true, "g")
% hold on
% plot(CL, CDwing, "r")
plot(CL, CDwhole6alt, "r")
xlabel("CL")
ylabel("CD")
grid on
yline(0)
xline(0)
legend("True", "Finite Wing", "Whole AC Approx", "Whole AC Approx (flat plate approx)")
title('Tempest C_D vs. C_L')

%% Functions
function out = altF(L, W, v, mu, rho) %summation of all diff shapes replaces CDmin
%inputs are length, width, velocity, and mu
Re_cr = 5*10^5;
xcr = Re_cr*mu/(rho*v);
q = 0.5*rho*v^2;
S = L*W;

if xcr > L
    Re_L = rho*v*L/mu;
    out = 1.328/sqrt(Re_L); %cf
else
    Re_Lwhole = rho*v*L/mu;
    Re_Lcrit = rho*v*xcr/mu;
    Cf_turb = .074/(Re_Lwhole^0.2)-.074/(Re_Lcrit^0.2);
    Cf_lamb = 1.328/sqrt(Re_Lcrit);
    Cf = Cf_turb+Cf_lamb;
    out = Cf;
end
end

function out = getvelo(CL,rho,W,S)
out = sqrt(W/(1/2*rho*S*CL));
end



    


