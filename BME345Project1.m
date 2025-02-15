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
diameterR1 = .102; % diameter of frame 1 on the bike (cm to m)
lengthCrank = .21; % m
height = .03; %m
width = .005; %m
%h1 = ; %height or length (m)
%d1 = ; %depth (m)
%h2 = ; %height or length (m)
%d2 = ; %depth (m)

% Mass of Links (kg)
m1 = 1.0; % frame mass approximation
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
figure(1)
subplot(3,1,1)
plot(th2new,th2new,"r",th2new,th3,"b",th2new,th4,"g","LineWidth",1.5)
title("\theta_2 \theta_3 \theta_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\theta_2, \theta_3, and \theta_4 in rads")
legend("\theta_2","\theta_3","\theta_4")

subplot(3,1,2)
plot(th2new,th3,"b")
plot(th2new,om2new,"r",th2new,om3,"b",th2new,om4,"g","LineWidth",1.5)
title("\omega_2 \omega_3 \omega_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\omega_2, \omega_3, and \omega_4 in rads/s")
legend("\omega_2","\omega_3","\omega_4")

subplot(3,1,3)
plot(th2new,al2new,"r",th2new,al3,"b",th2new,al4,"g","LineWidth",1.5)
title("\alpha_2 \alpha_3 \alpha_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\alpha_2, \alpha_3, and \alpha_4 in rads/(s^2)")
legend("\alpha_2","\alpha_3","\alpha_4")

%% Graph of Angle 2 vs. Torque of Pedal and Hips
figure(2)
plot(th2new,th2new,"r",th2new,th3,"b",th2new,th4,"g","LineWidth",1.5)
title("\theta_2 \theta_3 \theta_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\theta_2, \theta_3, and \theta_4 in rads")
legend("\theta_2","\theta_3","\theta_4")

%% Graph of the Force Components vs. Angle 2
figure(3)
plot(th2new,th2new,"r",th2new,th3,"b",th2new,th4,"g","LineWidth",1.5)
title("\theta_2 \theta_3 \theta_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\theta_2, \theta_3, and \theta_4 in rads")
legend("\theta_2","\theta_3","\theta_4")


%% Video Generation
