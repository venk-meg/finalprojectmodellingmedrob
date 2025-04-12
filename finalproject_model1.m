%% Main Script
% Clear workspace and command window
clear;
clc;

%% Parameters
PLA_den = 1240;             % PLA density kg/m3
S = 0.05;                   % Length of each side link (m)
T = 0.085;                  % Length of the top link (m)
D = T;                      % Distance between pivots (m)
w_s = 0.005;                % Width of side link (m)
w_t = 0.0108;               % Width of top link (m)
m_s = PLA_den * w_s * S;    % Mass of each side link (kg)
m_t = PLA_den * w_t * T;    % Mass of the top link (kg)

% Spring geometry parameters (for stiffness calculation)
% (Adjust these parameters to modify the spring stiffness.)
p1 = 0.1;    % L spring1 (m)
q1 = 0.05;    % W spring1 (m)
r1 = 0.05;   % H spring1 (m)
p2 = 0.3;    % L spring2 (m)
q2 = 0.05;    % W spring2 (m)
r2 = 0.05;   % H spring2 (m)
E = 55158;   % Elastic modulus (N/m)

% Compute spring stiffnesses
k1 = E * q1 * r1 / p1;   % Spring 1 stiffness (N/m)
k2 = E * q2 * r2 / p2;   % Spring 2 stiffness (N/m)

% For the equilibrium investigation, we will allow the spring rest lengths to vary.
% (Other parameters like attachment points could also be swept.)
l1_o = 0.1;             % initial rest length of spring 1 (m)
l2_o = 0.3;             % initial rest length of spring 2 (m)

% Attachment points along the side links
a = 0.02;                % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.02;                % Distance from bottom-right pivot to spring 2 attachment (m)

% Fixed bottom attachment points for springs
Q1 = -0.03;              % Spring 1 bottom attachment point (m)
Q2 = D + 0.015;          % Spring 2 bottom attachment point (m)

% Define the applied torque for the dynamic simulation; here set to zero.
tau = @(t) (-40).*(t < 1);
%tau = @(t) 0;

% Damping coefficient (added to reduce oscillations)
c = 2;   % Adjust this value as needed

% Initial conditions (start from rest)
theta0 = deg2rad(0);           % Initial angle (rad)
theta_dot0 = 0;          % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for the simulation
t_span = [0, 10];

%% Modular Functions for Spring Lengths and Their Derivatives
l1_fun = @(theta) sqrt((a * sin(theta))^2 + (a * cos(theta) + Q1)^2);
dl1_fun = @(theta) (a * cos(theta) + Q1) / l1_fun(theta);

l2_fun = @(theta) sqrt((b * sin(theta))^2 + ((Q2 - T) + b * cos(theta))^2);
dl2_fun = @(theta) ((Q2 - T) + b * cos(theta)) / l2_fun(theta);

%% Solve the ODE using ode45 (Dynamic Simulation)
[t, y] = ode45(@(t,y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, c, tau, l1_fun, dl1_fun, l2_fun, dl2_fun), t_span, initial_conditions);

%% Find Equilibrium Points using fzero
% Our goal is to find two equilibria (net torque = 0):
%   1) When theta is near 0 and applied (constant) torque tau_const = 0.
%   2) When theta is near 90° (pi/2 radians) and applied (constant) torque tau_const = -1.

tau_const1 = 0;
tau_const2 = -1;

theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const1), 0);
theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const2), pi/2);

disp(['Equilibrium for tau_const = 0 (target near theta = 0): ', num2str(theta_eq1), ' rad, ', num2str(rad2deg(theta_eq1)), ' deg']);
disp(['Equilibrium for tau_const = -1 (target near theta = 90 deg): ', num2str(theta_eq2), ' rad, ', num2str(rad2deg(theta_eq2)), ' deg']);

%% (Optional) Fancy Chart: Time-Series Plots
figure;
subplot(4,1,1);
plot(t, y(:,1), 'b-', 'LineWidth', 2);
ylabel('\theta (rad)');
title('Angle vs. Time');
grid on;

subplot(4,1,2);
plot(t, arrayfun(tau,t), 'm-', 'LineWidth', 2);
ylabel('Torque (N.m)');
title('Applied Torque vs. Time');
grid on;

% Precompute spring forces for plotting:
numPoints = length(t);
spring_force1 = zeros(numPoints,1);
spring_force2 = zeros(numPoints,1);
for i = 1:numPoints
    theta_i = y(i,1);
    spring_force1(i) = k1 * (l1_fun(theta_i) - l1_o) * dl1_fun(theta_i);
    spring_force2(i) = k2 * (l2_fun(theta_i) - l2_o) * dl2_fun(theta_i);
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

%% Animation with Live Data Display (optional)
figure;
hold on;
axis equal;
xlim([-0.1, 0.15]);
ylim([0, 0.2]);
xlabel('X (m)');
ylabel('Y (m)');
title('Mechanism Animation with Live Data');
plot([0, D], [0, 0], 'k-', 'LineWidth', 2);  % Fixed base

h_side_left = plot([0, 0], [0, 0], 'b-', 'LineWidth', 3);   % Left side link
h_side_right = plot([D, D], [0, 0], 'r-', 'LineWidth', 3);    % Right side link
h_top = plot([0, D], [0, 0], 'g-', 'LineWidth', 3);           % Top link

plot(0, 0, 'ko', 'MarkerFaceColor', 'k');  % Pivot markers
plot(D, 0, 'ko', 'MarkerFaceColor', 'k');

h_spring1 = plot([0,0], [0,0], 'b--');  % Spring 1 line
h_spring2 = plot([D,D], [0,0], 'r--');  % Spring 2 line
h_data = text(0.95, 0.9, '', 'FontSize', 12, 'Color', 'm', 'Units', 'normalized','HorizontalAlignment','right');

for i = 1:length(t)
    theta1 = y(i,1);
    
    % Compute free endpoint positions of side links (parallelogram constraint)
    A_x = -S * cos(theta1);
    A_y = S * sin(theta1);
    B_x = D - S * cos(theta1);
    B_y = S * sin(theta1);
    
    % Compute spring attachment points along side links:
    a_x = -a * cos(theta1);
    a_y = a * sin(theta1);
    b_x = D - b * cos(theta1);
    b_y = b * sin(theta1);
    
    set(h_side_left, 'XData', [0, A_x], 'YData', [0, A_y]);
    set(h_side_right, 'XData', [D, B_x], 'YData', [0, B_y]);
    set(h_top, 'XData', [A_x, B_x], 'YData', [A_y, B_y]);
    
    set(h_spring1, 'XData', [Q1, a_x], 'YData', [0, a_y]);
    set(h_spring2, 'XData', [Q2, b_x], 'YData', [0, b_y]);
    
    data_str = sprintf('Theta: %.1f deg\nTorque: %.2f N.m\nSpring1: %.2f N\nSpring2: %.2f N', ...
        rad2deg(theta1), tau(t(i)), spring_force1(i), spring_force2(i));
    set(h_data, 'String', data_str);
    
    pause(0.1);
end

%% Local Function Definitions
function dydt = system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, c, tau, l1_fun, dl1_fun, l2_fun, dl2_fun)
    theta1 = y(1);
    theta_dot1 = y(2);
    
    % Endpoint positions (parallelogram constraint: theta2 = theta1)
    A_x = -S * cos(theta1);
    A_y = S * sin(theta1);
    B_x = D - S * cos(theta1);
    B_y = S * sin(theta1);
    
    % Compute spring lengths and their derivatives via modular functions
    l1 = l1_fun(theta1);
    l2 = l2_fun(theta1);
    dl1 = dl1_fun(theta1);
    dl2 = dl2_fun(theta1);
    
    % Compute spring forces: F = k*(l - l0)*(dl/dθ)
    spring_force_1 = k1 * (l1 - l1_o) * dl1;
    spring_force_2 = k2 * (l2 - l2_o) * dl2;
    
    % External (applied) torque as function of time
    torque = tau(t);
    
    % Effective inertia (simplified lumped model)
    I_eff = (2/3)*m_s*S^2 + m_t*S^2;
    
    % Equation of motion: I_eff * theta_ddot = (spring forces) - (external torque)
    theta_ddot = (spring_force_1 + spring_force_2 - torque - c*theta_dot1) / I_eff;
    dydt = [theta_dot1; theta_ddot];
end

function netTorque = equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const)
    % Compute spring forces at a given theta
    l1 = l1_fun(theta);
    l2 = l2_fun(theta);
    dl1 = dl1_fun(theta);
    dl2 = dl2_fun(theta);
    
    springForce1 = k1 * (l1 - l1_o) * dl1;
    springForce2 = k2 * (l2 - l2_o) * dl2;
    
    % Equilibrium is reached when the net torque equals zero:
    % tau_const - (springForce1 + springForce2) = 0
    netTorque = tau_const - (springForce1 + springForce2);
end
