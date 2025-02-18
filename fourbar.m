function errors = fourbar(param,r1,r2,r3,r4,th1,th2,om2,al2)
%
% This function determines the error in the angular position, velocity, and
% acceleration equations for a four-bar linkage with a fixed, horizontal
% frame at r1. It is intended to be used in conjuction with fsolve to
% numerically solve for positions, velocities, and accelerations of
% non-input linkages.
%
% Input:
%     param: Initial guesses for angular position, velocity, and 
%            accelerations of linkages r3 and r4. Angles in radians.
%     r1, r2, r3, r4: Lengths of the linkages.
%     th2, om2, al2: Input parameters from link 2. Angle in radians
%
% Output:
%     errors: Errors in position, velocity, and accelerations. If the
%             inputs satisfy the three sets of equations, then the errors
%             vector would be zeros.

th3 = param(1);
th4 = param(2);
om3 = param(3);
om4 = param(4);
al3 = param(5);
al4 = param(6);

% Position Equations:
px = r1*cos(th1)+r2*cos(th2)+r3*cos(th3)+r4*cos(th4);
py = r1*sin(th1)+r2*sin(th2)+r3*sin(th3)+r4*sin(th4);

% Velocity Equations
vx = -om2*r2*sin(th2)-om3*r3*sin(th3)-om4*r4*sin(th4);
vy = om2*r2*cos(th2)+om3*r3*cos(th3)+om4*r4*cos(th4);

% Acceleration Equations
ax = (-al2*r2*sin(th2)-(om2^2)*r2*cos(th2))+(-al3*r3*sin(th3)-(om3^2)*r3*cos(th3))...
+ (-al4*r4*sin(th4)-(om4^2)*r4*cos(th4));
ay = (al2*r2*cos(th2)-(om2^2)*r2*sin(th2)) +(al3*r3*cos(th3)-(om3^2)*r3*sin(th3))...
    + (al4*r4*cos(th4)-(om4^2)*r4*sin(th4));

errors = [px py vx vy ax ay];