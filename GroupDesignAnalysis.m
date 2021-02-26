%% ----------------------- Team Glider Design -------------------------- %%
%% Cleaning
clear;
clc;
close all;
%% Variables
% Provided flat plate CL and CD values with respect to AOA
plateCl = [-.75 -.76 -.77 -.78 -0.8 -0.78 -.76 -.74 -.7 -.6 -.5 -.38 -.3 -.2 -.1 0 .1 .2 .3 .4 .49 .6 .7 .75 .76 .77 .76 .75 .75 .74 .74 .73 .72 .71 .71 .7];
plateCd = [.21 .2 .19 .18 .17 .15 .13 .11 .09 .07 .05 .04 .02 .01 .005 0  .005 .01 .015 .02 .05 .065 0.09 .11 .13 .15 .17 .19 .199 .21 .22 .23 .25 .27 .285 .29];
plateAOA = [-15:1:20];
ReCr = 140000;

% Air properties assuming altiitude of 1500 m
mu      = 1.74 * 10^-5;                                                    % Air viscosity [kg/m*s]
rhoAlt  = 1.0581;                                                          % Air density [kg/m^3]
ReL     = 1000000;                                                         % Local Reynold's value (can't remember where this number is from - maybe milestone 1
rhoFoam = 0.295 / 1000;                                                    % Density of Dollar Tree foam [kg/m^3]
g       = 9.81;                                                            % Acceleration of gravity [m/s^2]

% Design properties
h        = 7;                                                              % Initial launch height[m]
b        = .98;                                                            % Wing span [m]
c        = nan;                                                            % Chord of wing [m]
AR       = 6;                                                            % Aspect ratio
Swet     = nan;                                                            % Wetted area [m^2]
Sref     = 2 * c * b;                                                            % Wing planform area [m^2]
Srest    = nan;                                                            % Surface area of fuselage and tail [m^2]
Sh       = nan;                                                            % Planform area of horizontal tail [m^2]
cBar     = Sref / b;                                                       % Mean aerodynamic chord [m]
cr       = nan;                                                            % Root chord [m]
ct       = nan;                                                            % Tip chord [m]
lambda   = ct / cr;                                                        % Taper ratio
mass     = nan;                                                            % Weight [kg]
e        = 0.9;                                                            % Estimated span efficiency factor
% e0       = (1.78*(1-(0.045*(AR^0.68)))) - 0.64;                            % Oswald's efficiency factor
e0       = 1;                                                              % Oswald's efficiency factor
V        = 7;                                                              % Velocity [m/s]
Wfuse    = nan;                                                            % Fuselage weight [N]
WHTail   = nan;                                                            % Horizontal tail (?) weight [N]
WVTail   = nan;                                                            % Vertical tail (?) weight [N]
WBallast = nan;                                                            % Ballast weight [N]
WPay     = nan;                                                            % Payload weight [N]
Wwing    = Sref*rhoFoam*g;                                                 % Wing weight [N]
Wrest    = Wfuse + WHTail + WVTail + WBallast + WPay;                      % Weight of glider without wings [N]
% WTO      = Wfuse + Wwing + WHTail + WVTail + WBallast + WPay;              % Takeoff weight [N]
WTO      = Wrest + Wwing;                                                  % Takeoff weight [N]
WingLoad = WTO / Sref;                                                     % Wing loading [N]
EPitCon  = nan;                                                            % Elevator pitch control [deg] (+/- 10 deg)
xCg      = nan;                                                            % Longitudinal stability center of gravity location [m] (0c < xCg < 1.0c)
xACv     = nan;                                                            % Vertical tail aerodynamic center [m]
Vh       = nan;                                                            % Horizontal tail volume [m^3] (0.3 < Vh < 0.6)
Vv       = nan;                                                            % Vertial tail volume [m^3] (0.02 < Vv < 0.05)
B        = nan;                                                            % Spiral parameter (B > 5)
dAngle   = nan;                                                            % Dihedral angle [deg]
ep0      = 0;                                                              % Downwash on tail at L = 0
iT       = nan;                                                            % Horizontal tail incident angle [deg]
angTail  = nan;                                                            % Tail angle of attack [deg]
Cmcg     = nan;                                                            % Coefficient of moment about center of gravity (=0 for trimmed flight)
Cmacw    = nan;                                                            % Coefficient of moment about wing aerodynamic center (approximate with 2D airfoil value)
CLw      = nan;                                                            % Coefficient of lift for wing at max range or max endurance (= to whole aircraft CL)
CLht     = nan;                                                            % Coefficient of lift for horizontal tail (flat plate approx)

%% ----------------------------- Analysis ------------------------------ %%
%% Estimated Lift Curve
fit2D = polyfit(plateAOA, plateCl, 1);                                     % Best fit line for Cl vs AOA
a0 = fit2D(1);                                                             % Save value of slope
fit2D = polyfit(plateCl, plateAOA, 1);                                     % Best fit line for AOA vs Cl (for finding AOA @ L=0)
angle0lift = polyval(fit2D, 0);                                            % Angle of attack at 0 lift
a = a0 / (1 + (57.3 * a0)/(pi * e * AR));                                  % 3D lift slope
calcCL = a .* (plateAOA - angle0lift);                                     % 3D coefficient of lift

%% Whole Aircraft Drag Polar
xCr = (ReCr * mu) / (rhoAlt * V);
CfTurbTotal = 0.074/(ReL^0.2);                                             % Total Wing Turbulent Cf
CfTurbLeft = 0.0592/(ReCr^0.2);                                            % Turbulent Cf to be subtracted from total
CfLamLeft = 0.664/sqrt(ReCr);                                              % Laminar Cf
CfTotal = 2*(CfTurbTotal - CfTurbLeft + CfLamLeft);                        % Total Cf for both wings

CDMin = CfTotal * (Swet / Sref);
CDwing = plateCd + (calcCL.^2 / (pi * e * AR));
fitDrag = polyfit(CDwing, calcCL,2 );                                      % Best fit line for finding CLminD
CLminD = polyval(fitDrag, min(CDwing));                                    % CL at min drag
CD0 = CDMin + (CLminD^2 / (pi*e*AR));                                      % Parasite drag
CD = CD0 + (calcCL.^2 / (pi*e*AR));                                        % Calculate whole aircraft coefficient of drag

%% Weight/Mass Estimation and Sizing


%% Longitudinal Static Stability and Trim


%% Lateral Directional Stability


%% Estimated Performance
LoD = calcCL ./ CD;                                                        % Find L/D values
LoDMax = max(LoD);                                                         % Find maximum L/D value
maxCL = calcCL(LoD == LoDMax);                                             % CL at L/D max
maxCD = CD(LoD == LoDMax);                                                 % CD at L/D max
theta = maxCD / maxCL;                                                     % Glide angle [deg]
sinkRate = V * sin(theta);                                                 % Sink rate of glider [m/s]

maxRange = h * LoDMax;                                                     % Max glide range [m]
maxRVel = sqrt(2*WingLoad ./ (maxCL * rhoAlt * Sref));                     % Max range velocity [m/s]
maxREndur = h / sinkRate;                                                  % Glide endurance [s]

%% Graphs
% Lift Curve Plots
figure();
plot(plateAOA, plateCl)
hold on;
plot(plateAOA, calcCL)
yline(0);
xline(0);
hold off;
grid on;
title('Estimated Lift Curves')
xlabel('Angle of Attack [deg]')
ylabel('Coefficient of Lift')
legend('2D Approx.', '3D Approx.', 'Location', 'Northwest')

% Whole Aircraft Drag Polar
figure();
plot(plateAOA, CD)
grid on;
title('Whole Aircraft Drag Polar')
xlabel('Angle of Attack [deg]')
ylabel('Coefficient of Drag')

%% Outputs
cost = mass * 1000;
summary = table(maxRange, maxRVel, maxREndur, EPitCon, xCg, Vh, Vv, B, b, cost);
summary.Properties.VariableNames = {'MaxGlideRange', 'MaxGlideRangeVelocity', 'MaxGlideEndurance', 'ElevatorPitchControl', 'XcgLocation', 'HorzTailVolume', 'VertTailVolume', 'Spiral', 'Wingspan', 'Cost'};