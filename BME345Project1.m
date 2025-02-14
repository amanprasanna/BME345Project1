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

% Total height
%h = ; %(m)

%Gravity
g = -9.81; %m/s^2

% Lengths (m)
r1 = 0.665; %bike pedal to hip
r2 = 0.185; %foot to bike pedal
r3 = 0.44; %knee to foot 
r4 = 0.51; %hip to knee


% Link information (might not need these)
p = 985; %density (kg/m^3)
%l = ; %height or length (m)
%d = ; %depth (m)

% Mass Ratios from Anthropometric Table

% Length Ratios from Anthropometric Table

% Mass of Links (kg)
m1 = ;
m2 = ;
m3 = 0.061 * m; %foot and leg (knee to foot) using anthropometric table
m4 = 0.100 * m; %thigh (hip to knee) using anthropometric table

% Moments of Inertia (kg*m^2)
I1 = ;
I2 = ;
I3 = ;
I4 = ;

% Crank Inputs if any
Ta = ; %N*m
th2 =  * pi/180; %rads
w2 = ; %rad/s
al2 = ; %rad/s^2
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

%% Plotting



%% Video Generation
