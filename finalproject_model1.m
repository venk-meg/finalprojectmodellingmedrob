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

%{
% Attachment points along the side links
a = 0.025;     % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.05;     % Distance from bottom-right pivot to spring 2 attachment (m)

% Fixed bottom attachment points for springs
Q1 = -0.01;   % Spring 1 bottom attachment point (m)
Q2 = D + 0.01; % Spring 2 bottom attachment point (m)

% Spring geometry parameters (for stiffness calculation)

p1 = sqrt(a^2 + Q1^2)    % L spring1 (m)
q1 = 0.015;   % W spring1 (m)
r1 = 0.003;   % H spring1 (m)

p2 = sqrt(b^2 + (Q2-D)^2)  % L spring2 (m)
q2 = 0.01;   % W spring2 (m)
r2 = 0.002;   % H spring2 (m)

E = 55158;   % Elastic modulus (N/m)

% Compute spring stiffnesses
k1 = 2 * E * q1 * r1 / p1;   % Spring 1 stiffness (N/m)
k2 = 2 * E * q2 * r2 / p2;   % Spring 2 stiffness (N/m)

% Rest lengths for springs
l1_o = 0.02;    % Rest length of spring 1 (m)
l2_o = 0.04;    % Rest length of spring 2 (m)

%}

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

% External (applied) torque (function of time)
% theta0 = deg2rad(0); tau = @(t) -75 * ((exp(-t))^((exp(t))^1000)) -0.5;
% theta0 = deg2rad(90); tau = @(t) 75 * ((exp(-t))^((exp(t))^1000));
 theta0 = deg2rad(38); tau = @(t) 0;
% theta0 = deg2rad(90); tau = @(t) -0.5;

% Damping coefficient
c = 0.3;   % Adjust as needed

% Initial conditions: [theta, theta_dot]
%theta0 = deg2rad(70);  % Initial angle (rad)
theta_dot0 = 0;        % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for dynamic simulation
t_span = [0, 5];

%% Solve the ODE using ode45 (Dynamic Simulation)
[t, y] = ode45(@(t,y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, c, tau), t_span, initial_conditions);

%% Find Equilibrium Points using fzero
% Two cases:
% 1) Target equilibrium near theta = 0 with tau_const = 0.
% 2) Target equilibrium near theta = 90° (pi/2 rad) with tau_const = -1.
tau_const1 = 0;
tau_const2 = -0.5;

theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, tau_const1), 0);
theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, tau_const2), pi/2);

disp(['Equilibrium for tau_const = 0 (target near theta = 0): ', num2str(theta_eq1), ' rad, ', num2str(rad2deg(theta_eq1)), ' deg']);
disp(['Equilibrium for tau_const = -0.5 (target near theta = 90 deg): ', num2str(theta_eq2), ' rad, ', num2str(rad2deg(theta_eq2)), ' deg']);

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
    % Compute net spring forces from individual contributions:
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

%% Animation with Live Data Display (optional, with video saving)
figure;
hold on;
axis equal;
xlim([-0.1, 0.15]);
ylim([-0.05, 0.15]);
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

% Create a VideoWriter object for saving the animation as an MP4
video_filename = 'mechanism_animation.mp4';
v = VideoWriter(video_filename, 'MPEG-4');  % 'MPEG-4' provides mp4 output
v.FrameRate = 20;   % Adjust frame rate if desired
open(v);

for i = 1:length(t)
    theta_i = y(i,1);
    
    % Compute free endpoint positions (parallelogram constraint)
    A_x = -S * cos(theta_i);
    A_y = S * sin(theta_i);
    B_x = D - S * cos(theta_i);
    B_y = S * sin(theta_i);
    
    % Compute spring attachment points along side links:
    a_x = -a * cos(theta_i);
    a_y = a * sin(theta_i);
    b_x = D - b * cos(theta_i);
    b_y = b * sin(theta_i);
    
    set(h_side_left, 'XData', [0, A_x], 'YData', [0, A_y]);
    set(h_side_right, 'XData', [D, B_x], 'YData', [0, B_y]);
    set(h_top, 'XData', [A_x, B_x], 'YData', [A_y, B_y]);
    
    set(h_spring1, 'XData', [Q1, a_x], 'YData', [0, a_y]);
    set(h_spring2, 'XData', [Q2, b_x], 'YData', [0, b_y]);
    
    data_str = sprintf('Theta: %.1f deg\nTorque: %.2f N.m\nSpring1: %.2f N\nSpring2: %.2f N', ...
        rad2deg(theta_i), tau(t(i)), spring_force1(i), spring_force2(i));
    set(h_data, 'String', data_str);
    
    drawnow;
    
    % Capture the frame and write it to the video file
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    pause(0.01);
end

close(v);  % Close the video file

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
    % Equilibrium when net torque (applied vs. spring forces) is zero
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
