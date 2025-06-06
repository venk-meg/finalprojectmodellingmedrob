%% Main Script
% Clear workspace and command window
clear;
clc;

%% Parameters
PLA_den = 1240;             % PLA density (kg/m^3)
S = 0.05;                   % Length of each side link (m)
T = 0.085;                  % Length of the top link (m)
D = T;                      % Distance between pivots (m)
w_s = 0.005;                % Width of side link (m)
w_t = 0.0108;               % Width of top link (m)
m_s = PLA_den * w_s * S;    % Mass of each side link (kg)
m_t = PLA_den * w_t * T;    % Mass of the top link (kg)

% Attachment points along the side links
a = 0.05;     % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.02;     % Distance from bottom-right pivot to spring 2 attachment (m)

% Fixed bottom attachment points for springs
Q1 = -0.02;   % Spring 1 bottom attachment point (m)
Q2 = D + 0.00; % Spring 2 bottom attachment point (m) 


% Rest lengths for springs
l1_o = 0.05;    % Rest length of spring 1 (m)
l2_o = 0.01;    % Rest length of spring 2 (m)

% Spring geometry parameters (for stiffness calculation)

p1 = sqrt(a^2 + Q1^2)    % L spring1 (m)
q1 = 0.015;   % W spring1 (m)
r1 = 0.004;   % H spring1 (m)

p2 = sqrt(b^2 + (Q2-D)^2)  % L spring2 (m)
q2 = 0.01;   % W spring2 (m)
r2 = 0.002;   % H spring2 (m)

E = 45000;   % Elastic modulus (N/m)

% Compute spring stiffnesses
k1 = 2 * E * q1 * r1 / p1;   % Spring 1 stiffness (N/m)
k2 = 2 * E * q2 * r2 / p2;   % Spring 2 stiffness (N/m)

%% Time-series Simulation (for reference)
% Here the external torque is defined as a constant.
theta0 = deg2rad(90); 
tau = @(t) -0.5;   % External torque
theta_dot0 = 0;    % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];
t_span = [0, 0.5];

[t, y] = ode45(@(t,y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, 0.2, tau), t_span, initial_conditions);

%% Find Equilibrium Points using fzero (for reference)
% Two cases:
% 1) Target equilibrium near theta = 0 with tau_const = 0.
% 2) Target equilibrium near theta = 90° (pi/2 rad) with tau_const = -0.5.
tau_const1 = 0;
tau_const2 = -0.5;

theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, tau_const1), 0);
theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, tau_const2), pi/2);

disp(['Equilibrium for tau_const = 0 (target near 25°): ', num2str(theta_eq1), ' rad, ', num2str(rad2deg(theta_eq1)), ' deg']);
disp(['Equilibrium for tau_const = -0.5 (target near 90°): ', num2str(theta_eq2), ' rad, ', num2str(rad2deg(theta_eq2)), ' deg']);

%% (Optional) Fancy Chart: Time-Series Plots
figure;
subplot(4,1,1);
plot(t, y(:,1), 'b-', 'LineWidth', 2);
ylabel('\theta (rad)');
title('Angle vs. Time');
grid on;

subplot(4,1,2);
plot(t, arrayfun(tau, t), 'm-', 'LineWidth', 2);
ylabel('Torque (N.m)');
title('Applied Torque vs. Time');
grid on;

% Precompute spring forces for plotting:
numPoints = length(t);
spring_force1 = zeros(numPoints,1);
spring_force2 = zeros(numPoints,1);
for i = 1:numPoints
    theta_i = y(i,1);
    [l1, dl1, l2, dl2] = springParameters(theta_i, a, b, Q1, Q2, T);
    spring_force1(i) = k1 * (l1 - l1_o) * dl1;
    spring_force2(i) = k2 * (l2 - l2_o) * dl2;
end

subplot(4,1,3);
plot(t, spring_force1, 'r--', 'LineWidth', 2);
ylabel('Spring Force 1 (N)');
title('Spring Force 1 vs. Time');
grid on;

subplot(4,1,4);
plot(t, spring_force2, 'g--', 'LineWidth', 2);
ylabel('Spring Force 2 (N)');
xlabel('Time (s)');
title('Spring Force 2 vs. Time');
grid on;

%% Bistability Analysis: Starting Angle vs. Equilibrium Angle (with τ = 0)
tau_bistable = @(t) 0;  % Set external torque to zero for bistability evaluation

% Define a range of starting angles (in degrees)
starting_angles_deg = linspace(-20, 110, 100);
final_equilibrium_deg = zeros(size(starting_angles_deg));

% Extend simulation time to allow the system to settle
t_span_bistability = [0, 2];

for i = 1:length(starting_angles_deg)
    theta0_local = deg2rad(starting_angles_deg(i));
    initial_conditions_local = [theta0_local; 0];
    [t_local, y_local] = ode45(@(t,y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, 0.2, tau_bistable), t_span_bistability, initial_conditions_local);
    theta_final = y_local(end, 1);
    final_equilibrium_deg(i) = rad2deg(theta_final);
end

figure;
plot(starting_angles_deg, final_equilibrium_deg, 'ko-', 'LineWidth', 2);
xlabel('Initial Angle (deg)');
ylabel('Final Equilibrium Angle (deg)');
title('Bistability Analysis: Initial vs. Equilibrium Angle with \tau = 0');
grid on;

%% Local Function Definitions

function dydt = system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, c, tau)
    theta = y(1);
    theta_dot = y(2);
    
    % Compute net torque from springs using helper function
    spring_net = springNetTorque(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T);
    
    % Effective inertia (lumped model)
    I_eff = (2/3)*m_s*S^2 + m_t*S^2;
    
    % Applied external torque
    ext_torque = tau(t);
    
    % Equation of motion: I_eff*theta_ddot = (spring_net - ext_torque - c*theta_dot)
    theta_ddot = (spring_net - ext_torque - c*theta_dot) / I_eff;
    dydt = [theta_dot; theta_ddot];
end

function netTorque = equilibrium_func(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, tau_const)
    % Compute net spring torque at a given angle
    spring_net = springNetTorque(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T);
    % Equilibrium when net spring torque (against constant external torque) is zero
    netTorque = tau_const - spring_net;
end

function [l1, dl1, l2, dl2] = springParameters(theta, a, b, Q1, Q2, T)
    % Compute spring 1 length and derivative
    l1 = sqrt((a*sin(theta))^2 + (a*cos(theta) + Q1)^2);
    dl1 = (a*cos(theta) + Q1) / l1;
    % Compute spring 2 length and derivative
    l2 = sqrt((b*sin(theta))^2 + ((Q2 - T) + b*cos(theta))^2);
    dl2 = ((Q2 - T) + b*cos(theta)) / l2;
end

function net = springNetTorque(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T)
    % Compute the net torque produced by the two springs
    [l1, dl1, l2, dl2] = springParameters(theta, a, b, Q1, Q2, T);
    F1 = k1 * (l1 - l1_o) * dl1;
    F2 = k2 * (l2 - l2_o) * dl2;
    net = F1 + F2;
end
