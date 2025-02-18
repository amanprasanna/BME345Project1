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

% Distal Distances
thighDistal = 0.567; 
thighProximal = 0.433; 

lowerlegDistal = 0.394; 
lowerlegProximal = 0.606; 

% Link information (might not need these)
pHuman = 985; %density (kg/m^3)
pCarbonFiber = 1600; % density of carbon fiber bikes (kg/m^3)
diameterR1 = 0.102; % diameter of frame 1 on the bike (cm to m)
lengthCrank = 0.21; % m
heightCrank = .03; %m
width = .005; %m
%h1 = ; %height or length (m)
%d1 = ; %depth (m)
%h2 = ; %height or length (m)
%d2 = ; %depth (m)

% Mass of Links (kg)
m1 = 1.0; % frame mass approximation
m2 = .168 + .698; % pedal + pedal crank (kg)
m3 = 0.061*m; %foot and leg (knee to foot) using anthropometric table
m4 = 0.100*m; %thigh (hip to knee) using anthropometric table

% Radius of Gyration (w.r.t CoG) of Links
k1 = diameterR1 / sqrt(2); %hip to bike pedal
k2 = heightCrank / sqrt(12); %bike pedal to foot (assume center to bike pedal is a rectangular prism...need to take measurements on this)
k3 = 0.416; %foot and leg (from anthropometric table)
k4 = 0.323; %thigh (from anthropometric table)

% Moments of Inertia (kg*m^2)
I2y = m*(lengthCrank^2 + heightCrank^2) / 12; % rectangular prism
I3 = m3*k3^2;
I4 = m4*k4^2;

% Crank Inputs if any
Ta = 125; %N*m
th2 = 180*pi/180; %rads
om2 = 1; %rad/s
al2 = 0; %rad/s^2
th1 = 225*pi/180; %rads

% Forces Due to Gravity
F1g = m1*g; %N
F2g = m2*g; %N
F3g = m3*g; %N
F4g = m4*g; %N


% Steps 
stepSize = pi/12; 
maxRev = 3*pi; 
th2new = th2:stepSize:maxRev; 

integerArray = numel(th2new);

% T hip
Thip = ones(1,length(th2new)); 
[~, indexT0] = min(abs(th2new - (3*pi/2))); 
[~, indexT1] = min(abs(th2new - (5*pi/2))); 

Thip = Thip*Ta; 
Thip(indexT0:indexT1) = 0; 

% Guess in the form [th3 th4 om3 om4 al3 al4] with radians, not degrees
guess = [(pi/4) (pi/4) 1 1 1 1];
options = optimoptions('fsolve','Display','final'); 

%% Calculations

% Finding the th3, th4, om3, om4, al3, al4 based on initial and new guesses
for k = 1:integerArray
    ans = fsolve(@fourbar_rev1,guess,options,r1,r2,r3,r4,th1,th2new(k),om2,al2);
    %
    % Redefine the guess as the recently solved parameters
    guess = [ans(1) ans(2) ans(3) ans(4) ans(5) ans(6)];
    %
    % Record newly calculated parameters in a vector or matrix.
    th3(k) = ans(1);
    th4(k) = ans(2);
    om3(k) = ans(3);
    om4(k) = ans(4);
    al3(k) = ans(5);
    al4(k) = ans(6);
    
end

% Finding the x and y positions and al components based on th2, th3, and th4
for k = 1:length(th2new)
    r12x(k) = -(r2/2).*cos(th2new(k)); 
    r12y(k) = -(r2/2).*sin(th2new(k)); 
    r32x(k) = (r2/2).*cos(th2new(k)); 
    r32y(k) = (r2/2).*sin(th2new(k)); 
    
    r23x(k) = -(r3*lowerlegDistal).*cos(th3(k)); % lower leg cannot be divided by 2
    r23y(k) = -(r3*lowerlegDistal).*sin(th3(k)); 
    r43x(k) = (r3*lowerlegProximal).*cos(th3(k)); 
    r43y(k) = (r3*lowerlegProximal).*sin(th3(k)); 

    r34x(k) = -(r4*thighDistal).*cos(th4(k)); % thigh cannot be divided by 2 
    r34y(k) = -(r4*thighDistal).*sin(th4(k));
    r14x(k) = (r4*thighProximal).*cos(th4(k)); 
    r14y(k) = (r4*thighProximal).*sin(th4(k));

    al2x(k) = 0 - (om2).*(r2/2).*cos(th2new(k)); 
    al2y(k) = 0 - (om2).*(r2/2).*sin(th2new(k)); 
    
    al3x(k) = (0 - (om2).*(r2).*cos(th2new(k))) + (-al3(k).*(r3*lowerlegDistal).*sin(th3(k)) - (om3(k)^2)*(r3*lowerlegDistal)*cos(th3(k))); 
    al3y(k) = (0 - (om2).*(r2).*sin(th2new(k))) + (al3(k).*(r3*lowerlegDistal).*cos(th3(k)) - (om3(k)^2).*(r3*lowerlegDistal).*sin(th3(k)));
    
    al4x(k) = (0 - (om2).*(r2).*cos(th2new(k))) + (-al3(k).*(r3).*sin(th3(k)) - (om3(k)^2)*(r3)*cos(th3(k))) ...
    + (-al4(k).*(r4*thighDistal).*sin(th4(k)) - (om4(k).^2)*(r4*thighDistal).*cos(th4(k))); 
    al4y(k) = (0 - (om2).*(r2).*sin(th2new(k))) + (al3(k).*(r3).*cos(th3(k)) - (om3(k)^2).*(r3).*sin(th3(k))) ...
    + (al4(k).*(r4*thighDistal).*cos(th4(k)) - (om4(k).^2).*(r4*thighDistal).*sin(th4(k))); 

end


% Solving for the Force Components

%  F12x  F12y  F32x F32y F23x F23y F43x F43y F34x  F34y F14x F14y  TD
for k = 1:length(al4y)
A = [1    0     1    0    0    0    0     0    0    0    0    0    0; 
     0    1     0    1    0    0    0     0    0    0    0    0    0;
 -r12y(k)' r12x(k)' -r32y(k)' r32x(k)'  0     0    0    0    0    0    0    0    1; % torque to solve is T2
     0    0     0    0    1    0    1     0    0    0    0    0    0;
     0    0     0    0    0    1    0     1    0    0    0    0    0;
     0    0     0    0  -r23y(k)' r23x(k)' -r43y(k)' r43x(k)'  0    0    0    0    0;
     0    0     0    0    0    0    0     0    1    0    1    0    0;
     0    0     0    0    0    0    0     0    0    1    0    1    0;
     0    0     0    0    0    0    0     0  -r34y(k)' r34x(k)' -r14y(k)' r14x(k)' 0;
     0    0     1    0    1    0    0     0    0    0    0    0    0;
     0    0     0    1    0    1    0     0    0    0    0    0    0;
     0    0     0    0    0    0    1     0    1    0    0    0    0;
     0    0     0    0    0    0    0     1    0    1    0    0    0];

% % %Right hand side of the equation/ known variables
b = [m2.*al2x(k) (m2.*al2y(k) - m2*g) (I2y*al2) m3.*al3x(k) (m3.*al3y(k) - m3*g) I3.*al3(k) m4.*al4x(k) (m4.*al4y(k) - m4*g) (I4.*al4(k) - Thip(k)) 0 0 0 0]';

% % % solving the system of Equations
F(:,k)= A\b;

end

% Solving for Thigh Force Components to Graph 
FthighParallel = F(end-5, :).*cos(th4); % force of the thigh in x direction (N)
FthighPerpendicular = F(end-4, :).*cos(th4 + (pi/2)); % force of the thigh in y direction (N)

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
plot(th2new,om2,"r",th2new,om3,"b",th2new,om4,"g","LineWidth",1.5)
title("\omega_2 \omega_3 \omega_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\omega_2, \omega_3, and \omega_4 in rads/s")
legend("\omega_2","\omega_3","\omega_4")

subplot(3,1,3)
plot(th2new,al2,"r",th2new,al3,"b",th2new,al4,"g","LineWidth",1.5)
title("\alpha_2 \alpha_3 \alpha_4 vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("\alpha_2, \alpha_3, and \alpha_4 in rads/(s^2)")
legend("\alpha_2","\alpha_3","\alpha_4")

%% Graph of Angle 2 vs. Torque of Pedal and Hips
figure(2)
plot(th2new,Thip,"r",th2new,F(end, :),"b","LineWidth",1.5)
title("Torque at the Pedal and Hips vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("Torque (N*m)")
legend("Torque of Hips","Torque of Pedals")

%% Graph of the Force Components vs. Angle 2
figure(3)
plot(th2new,F4,"r",th2new,th2new,"b", "LineWidth",1.5)
title("Force Components Parallel & Perpendicular to the Thigh vs. \theta_2")
xlabel("\theta_2 in radians")
ylabel("Force (N)")
legend("Force Parallel to Thigh","Force Perpendicular to Thigh")


%% Video Generation
