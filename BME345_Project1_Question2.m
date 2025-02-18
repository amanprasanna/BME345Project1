% Four-bar linkage kinematics analysis

clear
clc
close all

%% Declarations
r1 = 0.665; r2 = 0.185; r3 = 0.44; r4 = 0.51;
om2 = 1;
al2 = 0;
th1 = (210*pi)/180;


% Define the range for theta2 (in radians)
revolutions = 2;
steps = 200; % Number of steps per revolution
th2_range = linspace(pi, 3*pi, steps);

% Initial guess for parameters of the third and fourth linkages
% [th3, th4, om3, om4, al3, al4]
guess = [(120*pi)/180, (20*pi)/180, 1, 1, 0, 0]; % Simplified initial guess

% Set options for fsolve
options = optimoptions('fsolve','Display','iter','MaxIter',1000,...
    'MaxFunctionEvaluations',5000,'FunctionTolerance',1e-10);

% Initialize arrays to store results
results = zeros(length(th2_range), 6);

%% Calculations
%% Calculations
for k = 1:length(th2_range)
    th2 = th2_range(k);  % Add semicolon
    
    % Calculate omega2 based on alpha2
    %om2 = sqrt(om2_initial^2 + 2*al2*(th2 - th2_initial));  % Add semicolon
    
    % Use fsolve to calculate new values
    ans(k,:) = fsolve(@(param) fourbar(param, r1, r2, r3, r4, th1, th2, om2, al2), guess, options);  % Add semicolon
    
    % Redefine the guess as the recently solved parameters
    guess = ans(k,:);  % Add semicolon
    
    % Record newly calculated parameters
    results = ans;  % Add semicolon
end



% Extract individual parameters from results
th3 = results(:, 1);
th4 = results(:, 2);
om3 = results(:, 3);
om4 = results(:, 4);
al3 = results(:, 5);
al4 = results(:, 6);

%% Plotting
figure(1)

% Theta plot
subplot(3,1,1)
plot(th2_range, th3, 'r-', ...
     th2_range, th4, 'g-')
title('Angular Position vs. \theta_2')
xlabel('\theta_2 (rad)')
ylabel('\theta (rad)')
legend('\theta_3', '\theta_4', 'Location', 'eastoutside')
grid on

% Omega plot
subplot(3,1,2)
plot(th2_range, om3, 'r-', ...
     th2_range, om4, 'g-')
title('Angular Velocity vs. \theta_2')
xlabel('\theta_2 (rad)')
ylabel('\omega (rad/s)')
legend('\omega_3', '\omega_4', 'Location', 'eastoutside')
grid on

% Alpha plot
subplot(3,1,3)
plot(th2_range, al3, 'r-', ...
     th2_range, al4, 'g-')
title('Angular Acceleration vs. \theta_2')
xlabel('\theta_2 (rad)')
ylabel('\alpha (rad/s^2)')
legend('\alpha_3', '\alpha_4', 'Location', 'eastoutside')
grid on

% Adjust subplot spacing
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Arial', 'FontSize', 10, 'XMinorTick', 'on', 'YMinorTick', 'on')


% Question 1 
rfx = 0.25*cosd(285);
rfy = 0.25*sind(285);
ruax = 0.28*cosd(155);
ruay = 0.28*sind(155);
ruaxcom = 0.564*0.28*cosd(155);
ruaycom = 0.564*0.28*sind(155);

alphaua = -1.5;
alphaf = 0.2;
omegaua = -8;
omegaf = 2;

vx = (-omegaf*rfy) + (-omegaua*ruay);
vy = (omegaf*rfx) + (omegaua*ruax);

Ax = (-alphaua*ruaycom - ((omegaua)^2)*ruaxcom)+(-alphaf*rfy - ((omegaf)^2)*rfx);

A = (alphaua*ruaxcom - ((omegaua)^2)*ruaycom) + (alphaf*rfx - ((omegaf)^2)*rfy);
