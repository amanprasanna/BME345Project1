% BME 345 Project 1
% Jasmine, John, Aman, Mylah
% Due: February 18, 2025

% Note: lots of stuff declared, may not need all of it

clear
clc
close all

%% Declarations

% Total mass
mHuman = 65.0; %(kg)

%Gravity
g = -9.81; %m/s^2

% Lengths (m)
r1 = 0.665; %hip to bike pedal
r2 = 0.185; %bike pedal to foot
r3 = 0.44; %foot to knee 
r4 = 0.51; %knee to hip


% Link information-Density (kg/m^3)
p1 = ; %density of carbon fiber
p2 = ; %density of carbon fiber
p3 = 985; %density of human body 
p4 = 985; %density of human body 

% Link information-Height/Length (m)
%h1 = ; 
%h2 = ;
%h3 = ;
%h4 = ;

% Link information-Depth (m)
%d1 = ; 
%d2 = ;
%d3 = ;
%d4 = ;


% Mass of Links (kg)
m1 = ;
m2 = ;
m3 = 0.061 * mHuman; %foot and leg (knee to foot) using anthropometric table
m4 = 0.100 * mHuman; %thigh (hip to knee) using anthropometric table

% Radius of Gyration (w.r.t CoG) of Links
k1 = ; %hip to bike pedal
k2 = ; %bike pedal to foot
k3 = 0.416; %foot and leg (from anthropometric table)
k4 = 0.323; %thigh (from anthropometric table)

% Moments of Inertia (kg*m^2)
I1 = m1 * k1^2;
I2 = m2 * k2^2;
I3 = m3 * k3^2;
I4 = m4 * k4^2;

% Crank Inputs
Ta4 = [0, 125]; %Applied Torque from Thigh in N*m

th2 = __ * pi/180; %rads
w2 = 1; %rad/s
al2 = 2; %rad/s^2

th1 =  * pi/180; %rads


% Forces Due to Gravity
F1g = m1 * g; %N
F2g = m2 * g; %N
F3g = m3 * g; %N
F4g = m4 * g; %N

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
