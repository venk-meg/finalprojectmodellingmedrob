% Clear workspace
clear;
clc;

% Parameters
PLA_den = 1240;         % PLA density kg/m3
S = 0.05;              % Length of each side link (m)
T = 0.085;             % Length of the top link (m)
w_s =  0.005;          % width side
w_t =  0.0102;         % width top
m_s = PLA_den * w_s * S;   % Mass of each side link (kg)
m_t = PLA_den * w_t * T;   % Mass of the top link (kg)
k1 = 10;               % Spring 1 stiffness (N/m)
k2 = 10;               % Spring 2 stiffness (N/m)
l1_o = 0.03;           % Rest length of spring 1 (m)
l2_o = 0.03;           % Rest length of spring 2 (m)
a = 0.025;             % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.025;             % Distance from bottom-right pivot to spring 2 attachment (m)
D = 0.085;             % Distance between pivots (m)
Q1 = -0.05;            % Spring 1 bottom attachment point (m)
Q2 = D + 0.05;         % Spring 2 bottom attachment point (m)

% Define the applied torque (smooth increase/decrease with theta)
tau = @(theta) 0.0 * sin(theta);  % Example applied torque (smooth)

% Initial conditions (start from rest)
theta0 = 0;             % Initial angle (rad)
theta_dot0 = 0;         % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for the simulation
t_span = [0, 20];       % Simulate from t=0 to t=20 seconds (slower simulation)

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
    
    % Spring forces (passive forces, only reacting to displacement)
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
axis equal;  % Keep aspect ratio fixed
xlim([-0.1, 0.15]);  % Adjust X-axis limits for better visualization
ylim([0, 0.2]);  % Adjust Y-axis limits for better visualization

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
h_spring1 = plot([0, 0], [0, 0], 'b--'); % Spring 1 -> between a and A
h_spring2 = plot([D, D], [0, 0], 'r--'); % Spring 2 -> between b and B

% Plot the torque indication as an arrow (optional)
h_torque_text = text(0.1, 0.18, '', 'FontSize', 12, 'Color', 'm');  % Placeholder text for torque

% Loop for animation
for i = 1:length(t)
    % Extract the angle and angular velocity from the solution
    theta1 = y(i, 1);
    
    % Calculate positions of key points
    A_x = -S * cos(theta1);  % Position of point A (left side link end)
    A_y = S * sin(theta1);
    
    B_x = D - S * cos(theta1);  % Position of point B (right side link end)
    B_y = S * sin(theta1);
    
    % Calculate spring attachment points a and b
    a_x = -a * cos(theta1);  % Position of spring attachment point a
    a_y = a * sin(theta1);  % Position of spring attachment point a
    
    b_x = D - b * cos(theta1);  % Position of spring attachment point b
    b_y = b * sin(theta1);      % Position of spring attachment point b
    
    % Update the positions of the links
    set(h_side_left, 'XData', [0, A_x], 'YData', [0, A_y]);
    set(h_side_right, 'XData', [D, B_x], 'YData', [0, B_y]);
    set(h_top, 'XData', [A_x, B_x], 'YData', [A_y, B_y]);
    
    % Update spring positions (for visualization)
    set(h_spring1, 'XData', [Q1, a_x], 'YData', [0, a_y]);  % Update Spring 1 -> between Q1 and a
    set(h_spring2, 'XData', [Q2, b_x], 'YData', [0, b_y]);  % Update Spring 2 -> between Q2 and b
    
    % Update torque arrow (optional)
    torque = tau(theta1);  % Get current torque value
    torque_magnitude = 0.03 * torque;  % Adjust the torque magnitude for visualization
    % Update torque value text
    torque = tau(theta1);  % Get current torque value
    set(h_torque_text, 'String', ['Torque: ' num2str(torque, '%.2f') ' N.m']);  % Update torque value displayed
    
    % Pause for animation effect (slower simulation)
    pause(0.1);  % Slower
end