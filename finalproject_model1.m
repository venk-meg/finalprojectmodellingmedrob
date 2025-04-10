% Clear workspace
clear;
clc;

% Parameters
PLA_den = 1240 ;         %PLA density kg/m3
S = 0.05;                   % Length of each side link (m)
T = 0.085;                   % Length of the top link (m)
w_s =  0.005                %width side
w_t =  0.0102               %width top
m_s = PLA_den*w_s*S;                 % Mass of each side link (kg) 0.00124g/mm3 for PLA * vol
m_t = PLA_den*w_t*T;                 % Mass of the top link (kg)
k1 = 10;                % Spring 1 stiffness (N/m)
k2 = 10;                % Spring 2 stiffness (N/m)
l1_o = 0.03;              % Rest length of spring 1 (m)
l2_o = 0.03;              % Rest length of spring 2 (m)
a = 0.025;                   % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.025;                   % Distance from bottom-right pivot to spring 2 attachment (m)
D = 2;                   % Distance between pivots (m)
Q1 = -0.05;                % Spring 1 bottom attachment point (m)
Q2 = D + 0.05;                % Spring 2 bottom attachment point (m)

% Define the applied torque (for example, sinusoidal)
%tau = @(theta) 10 * sin(theta);  % Example applied torque
tau = @(theta) 10 * 0.01;  % Example applied torque

% Initial conditions
theta0 = 0;              % Initial angle (rad)
theta_dot0 = 0;          % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for the simulation
t_span = [0, 10];        % Simulate from t=0 to t=10 seconds

% Differential equation function (returns derivatives)
function dydt = system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau)
    theta1 = y(1);  % Angle of the left side link
    theta_dot1 = y(2);
    
    % For the parallelogram constraint, we know theta2 = theta1
    % Calculate the positions of points A and B, with fixed length of the top link
    A_x = -S * cos(theta1);
    A_y = S * sin(theta1);
    
    B_x = D - S * cos(theta1); % Fixed top link length constraint
    B_y = S * sin(theta1);     % Symmetric behavior
    
    % Calculate spring lengths
    l1 = sqrt((a * sin(theta1))^2 + (a * cos(theta1) + Q1)^2); % Spring 1 length
    l2 = sqrt((b * sin(theta1))^2 + ((Q2 - T) + b * cos(theta1))^2); % Spring 2 length
    
    % Derivatives of spring lengths
    dl1_dtheta = (a * cos(theta1) + Q1) / l1;
    dl2_dtheta = ((Q2 - T) + b * cos(theta1)) / l2;
    
    % Spring forces
    spring_force_1 = k1 * (l1 - l1_o) * dl1_dtheta;
    spring_force_2 = k2 * (l2 - l2_o) * dl2_dtheta;
    
    % Applied torque
    torque = tau(theta1);
    
    % Effective inertia
    I_eff = (2 / 3) * m_s * S^2 + m_t * S^2 * cos(theta1)^2;
    
    % Differential equation for theta1 (since theta2 = theta1 in the parallelogram)
    theta1_ddot = (spring_force_1 + spring_force_2 - torque) / I_eff;
    
    % Return the derivatives
    dydt = [theta_dot1; theta1_ddot];
end

% Solve the ODE using ode45
[t, y] = ode45(@(t, y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau), t_span, initial_conditions);

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
    theta1 = y(i, 1);
    % Calculate positions of key points
    A_x = -S * cos(theta1);  % Position of point A (left side link end)
    A_y = S * sin(theta1);
    
    B_x = D - S * cos(theta1);  % Position of point B (right side link end)
    B_y = S * sin(theta1);
    
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
