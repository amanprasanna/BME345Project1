% BME 345 Project 1
% Jasmine, John, Aman, Mylah
% Due: February 18, 2025

% Note: lots of stuff declared, may not need all of it

clear
clc
close all

%% Declarations

% Total mass
m = 65.0; %(kg)

%Gravity
g = -9.81; %m/s^2

% Lengths (m)
r1 = 0.665; %hip to bike pedal
r2 = 0.185; %bike pedal to foot
r3 = 0.44; %foot to knee 
r4 = 0.51; %knee to hip


% Link information (might not need these)
pHuman = 985; %density (kg/m^3)
pCarbonFiber = 1600; % density of carbon fiber bikes (kg/m^3)
%diameterR1 = ; % diameter of frame 1 on the bike (cm to m)
%h1 = ; %height or length (m)
%d1 = ; %depth (m)
%h2 = ; %height or length (m)
%d2 = ; %depth (m)

% Mass of Links (kg)
m1 = ;
m2 = .168 + .698; % pedal + pedal crank (kg)
m3 = 0.061 * m; %foot and leg (knee to foot) using anthropometric table
m4 = 0.100 * m; %thigh (hip to knee) using anthropometric table

% Radius of Gyration (w.r.t CoG) of Links
k1 = diameterR1 / sqrt(2); %hip to bike pedal
k2 = ; %bike pedal to foot (assume center to bike pedal is a rectangular prism...need to take measurements on this)
k3 = 0.416; %foot and leg (from anthropometric table)
k4 = 0.323; %thigh (from anthropometric table)

% Moments of Inertia (kg*m^2)
I1 = m1 * k1^2;
I2 = m2 * k2^2;
I3 = m3 * k3^2;
I4 = m4 * k4^2;

% Crank Inputs if any
Ta = [0 120]; %N*m
th2 =  * pi/180; %rads
w2 = 1; %rad/s
al2 = 0; %rad/s^2
th1 =  * pi/180; %rads


% Forces Due to Gravity
F1g = ; %N
F2g = ; %N
F3g = ; %N
F4g = ; %N

% Guess in the form [th3 th4 om3 om4 al3 al4] with radians, not degrees
guess = [7*pi/4 3*pi/4 1 1 1 1];
options = optimoptions('fsolve','Display','final')

%% Calculations

ans = fsolve(@fourbar_rev1,guess,options,r1,r2,r3,r4,th1,th2,w2,al2)
th3 = ans(1);
th4 = ans(2); 
w3 = ans(3);
w4 = ans(4);
al3 = ans(5);
al4 = ans(6);

%% Graph of the Angles, Angular Accelerations, and Angular Velocities


%% Graph of Angle 2 vs. Torque of Pedal and Hips


%% Graph of the Force Components vs. Angle 2



%% Video Generation
