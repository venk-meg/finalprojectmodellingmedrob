% Clear workspace
clear;
clc;

% Parameters
m = 1;                   % Mass of each side link (kg)
m_top = 1;               % Mass of the top link (kg)
L = 1;                   % Length of each side link (m)
L_top = 2;               % Length of the top link (m)
k1 = 100;                % Spring 1 stiffness (N/m)
k2 = 100;                % Spring 2 stiffness (N/m)
l0_1 = 1.5;              % Rest length of spring 1 (m)
l0_2 = 1.5;              % Rest length of spring 2 (m)
D = 2;                   % Distance between pivots (m)
tau_act = @(theta) 10 * sin(theta); % Actuated torque (N.m) [example: sinusoidal torque]

% Define the system equations (ODE)
% State vector: [theta, theta_dot]
% theta is the generalized coordinate (angle)
% theta_dot is the velocity

% Initial conditions
theta0 = 0;              % Initial angle (rad)
theta_dot0 = 0;          % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for the simulation
t_span = [0, 10];        % Simulate from t=0 to t=10 seconds

% Define the system of ODEs
% ODEs are derived from the Lagrangian approach

% Differential equation function (returns derivatives)
function dydt = system_eqns(t, y, L, m, m_top, k1, k2, l0_1, l0_2, D, tau_act)
    theta = y(1);
    theta_dot = y(2);
    
    % Calculate spring lengths as functions of theta (modify as per your geometry)
    l1 = sqrt((L * cos(theta) - 0)^2 + (L * sin(theta))^2); % Spring 1 length
    l2 = sqrt((D - L * cos(theta))^2 + (L * sin(theta))^2); % Spring 2 length
    
    % Compute derivatives of spring lengths (using chain rule)
    dl1_dtheta = (L * cos(theta) - 0) / l1; % Derivative of spring length 1 with respect to theta
    dl2_dtheta = (D - L * cos(theta)) / l2; % Derivative of spring length 2 with respect to theta
    
    % Compute potential energy derivatives
    spring_force_1 = k1 * (l1 - l0_1) * dl1_dtheta; % Spring 1 generalized force
    spring_force_2 = k2 * (l2 - l0_2) * dl2_dtheta; % Spring 2 generalized force
    
    % Applied torque (as a function of theta)
    tau = tau_act(theta); % External applied torque
    
    % Equation of motion (Lagrangian formulation)
    % Assuming I is the inertia of side links and top link combined
    I_eff = (2 / 3) * m * L^2 + m_top * L^2 * cos(theta)^2;  % Effective inertia
    
    % Differential equations: [theta_dot, theta_ddot]
    theta_ddot = (spring_force_1 + spring_force_2 - tau) / I_eff;
    
    % Return the derivatives
    dydt = [theta_dot; theta_ddot];
end

% Solve the ODE using ode45
[t, y] = ode45(@(t, y) system_eqns(t, y, L, m, m_top, k1, k2, l0_1, l0_2, D, tau_act), t_span, initial_conditions);

% Prepare for animation
figure;
hold on;

% Plot the fixed base (bottom link)
plot([0, D], [0, 0], 'k-', 'LineWidth', 2);  % Base link

% Create plot for the mechanism animation
h_side_left = plot([0, 0], [0, 0], 'b-', 'LineWidth', 3);   % Left side link
h_side_right = plot([D, D], [0, 0], 'r-', 'LineWidth', 3);  % Right side link
h_top = plot([0, D], [0, 0], 'g-', 'LineWidth', 3);        % Top link

% Add markers for the joints (pivots)
h_left_pivot = plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Left pivot
h_right_pivot = plot(D, 0, 'ko', 'MarkerFaceColor', 'k'); % Right pivot

% Add spring lines (optional, for visualization)
h_spring1 = plot([0, 0], [0, 0], 'b--'); % Spring 1
h_spring2 = plot([D, D], [0, 0], 'r--'); % Spring 2

% Loop for animation
for i = 1:length(t)
    % Extract the angle and angular velocity from the solution
    theta = y(i, 1);
    % Calculate positions of key points
    A_x = L * cos(theta);  % Position of point A (left side link end)
    A_y = L * sin(theta);
    
    B_x = D - L * cos(theta);  % Position of point B (right side link end)
    B_y = L * sin(theta);
    
    % Update the positions of the links
    set(h_side_left, 'XData', [0, A_x], 'YData', [0, A_y]);
    set(h_side_right, 'XData', [D, B_x], 'YData', [0, B_y]);
    set(h_top, 'XData', [A_x, B_x], 'YData', [A_y, B_y]);
    
    % Update spring positions (optional for visualization)
    set(h_spring1, 'XData', [0, A_x], 'YData', [0, A_y]);
    set(h_spring2, 'XData', [D, B_x], 'YData', [0, B_y]);
    
    % Pause for animation effect
    pause(0.05);  % Adjust the pause for smoother or faster animation
end
