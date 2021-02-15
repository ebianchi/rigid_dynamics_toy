%%% Toy Rod Simulation
% Bibit Bianchini
% 2/14/2021

% This script aims to recreate the experiment in Stewart and Trinkle 1996
% in which a rigid rod with rounded ends falls onto a hard ground, and it
% experiences rotations, gravity, and friction.

clear;
clc;
close all;

%% Parameters
% System parameters
l = 0.5;   % [m]      - length of the rod, excluding the rounded ends
m = 1;     % [kg]     - mass of the rod
r = 0.05;  % [m]      - radius of the ends (half thickness of the rod)
I = 0.002; % [kg m^2] - moment of inertia
mu = 0.6;  %          - coefficient of friction between rod and table
g = -9.81; % [m/s^2]  - gravitational acceleration

% Simulation parameters
h = 0.0025;  % [s]    - step size
tf = 1;      % [s]    - final time for simulation

% Initial conditions
x0     = 0;            % [m]
y0     = 1;            % [m]
theta0 = deg2rad(30);  % [degrees]

dx0     = 0;  % [m/s]
dy0     = 0;  % [m/s]
dtheta0 = 4;  % [rad/s]

% Problem quantities
n = 3;   % number of generalized coordinates [x, y, theta]
p = 2;   % number of contacts
k = 2;   % friction polyhedron vertices (always 2 in 2D)


%% Static Vectors and Matrices
M = diag([m, m, I]);
Mu = diag([mu, mu]);
E = blkdiag(ones(k,1), ones(k,1));  % list ones(k,1) p times

% N, D, and k matrices/vectors are functions of q
% N     = N_rod(q_0, v_0, l);
% k_vec = k_rod(q_0, v_0, M);
% D     = D_rod(q_0, l);

%% Set up the configuration and speed matrices
ts = 0:h:tf;

qs = zeros(n, size(ts,2));
vs = zeros(n, size(ts,2));

qs(:,1) = [ x0;  y0;  theta0];
vs(:,1) = [dx0; dy0; dtheta0];

%% Iterate through time
for i = 1:size(ts,2)-1
    t = ts(i);
    q = qs(:,i);
    v = vs(:,i);
    
    % get the matrices/vectors that change as function of q, and
    % approximate q to be midpoint of next timestep
    k_vec = k_rod(q + h/2*v, v, M);
    N     = N_rod(q, v, l);
    D     = D_rod(q, l);

    phis = phis_rod(q, l, r);
    
    % build the large W matrix and b vector for the LCP form:
    % w = W*z + b,   w and z >= 0,   w'*z = 0
    W = [D' * (M \ D), D' * (M \ N), E;
         N' * (M \ D), N' * (M \ N), zeros(p,p);
         -E',          Mu,           zeros(p,p)];
    
    b = [D' * (v + h * (M\k_vec));
         phis/h + N' * (v + h * (M\k_vec));
         zeros(p,1)];
    % ^note:  there's a change in the second line to swap the less
    % intuitive constraint N' * q >= alpha from the Stewart and Trinkle
    % paper to the contactnets' constraint phrasing phi >= 0.
     
    % solve the LCP
    [w,z,retcode] = LCPSolve(W,b);
    
    % extract the Beta and Cn vectors from the z vector
    Beta = z(1:p*k);
    Cn = z(p*k+1 : p*k+p);
    
    % solve for the next states
    vs(:,i+1) = vs(:,i) + M \ (N*Cn + D*Beta + h*k_vec);
    qs(:,i+1) = qs(:,i) + h*vs(:,i+1);
end

%% Visualize the system dynamics
% Given the evolution of the generalized coordinates (qs) and the
% generalized speeds (vs) over time (ts), generate some helpful plots.

% Plot the velocity profile (Fig. 3)
figure();
plot(ts, vs(1,:));
hold on;
plot(ts, vs(2,:));
plot(ts, vs(3,:));
xlabel('Time (s)');
ylabel('Velocity (m/s or rad/s)');
legend('X Velocity', 'Y Velocity', '\theta Velocity');
title('Generalized Coordinate Velocities');

% Plot the positions
figure();
plot(ts, qs(1,:));
hold on;
plot(ts, qs(2,:));
plot(ts, qs(3,:));
xlabel('Time (s)');
ylabel('Generalized Coordinates (m or rad)');
legend('X', 'Y', '\theta');
title('Generalized Coordinates');

% Draw the central line of the rod over time
lefts  = [qs(1,:) - l/2 * cos(qs(3,:));
          qs(2,:) - l/2 * sin(qs(3,:))];
rights = [qs(1,:) + l/2 * cos(qs(3,:));
          qs(2,:) + l/2 * sin(qs(3,:))];
figure();
hold on;
for i = 1:size(lefts,2)
    plot([lefts(1,i), rights(1,i)], [lefts(2,i), rights(2,i)]);
end
axis equal;
title('Rod Over Time');


