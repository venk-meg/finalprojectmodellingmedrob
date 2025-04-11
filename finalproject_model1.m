% Clear workspace
clear;
clc;

% Parameters
PLA_den = 1240;         % PLA density kg/m3
S = 0.05;              % Length of each side link (m)
T = 0.085;             % Length of the top link (m)
D = T;                 % Distance between pivots (m)
w_s =  0.005;          % width side
w_t =  0.0108;         % width top
m_s = PLA_den * w_s * S;   % Mass of each side link (kg)
m_t = PLA_den * w_t * T;   % Mass of the top link (kg)

p1 = 0.4 ;    %L spring1 ,m
q1 = 0.2  ; %W spring1 ,m
r1 =0.05   ;     %H spring1 ,m
p2  = 0.4   ;    %L spring2 ,m
q2  = 0.2    ;   %W spring2 ,m
r2   =0.05    ;   %H spring2 ,m

E =  55158      ; %Elastic modulus estimate (linear for now nonlinear later), 55158 N/m which is 8 psi

k1 = E*q1*r1/p1;               % Spring 1 stiffness (N/m)
k2 = E*q2*r2/p2;               % Spring 2 stiffness (N/m)
l1_o = 0.07;           % Rest length of spring 1 (m)
l2_o = 0.03;           % Rest length of spring 2 (m)
a = 0.05;             % Distance from bottom-left pivot to spring 1 attachment (m)
b = 0.05;             % Distance from bottom-right pivot to spring 2 attachment (m)
Q1 = -0.03;            % Spring 1 bottom attachment point (m)
Q2 = D + 0.015;         % Spring 2 bottom attachment point (m


% Define the applied torque (smooth increase from 0 to 1 over time)
%tau = @(t) ((mod(t,4) >= 1 & mod(t,4) < 2) * (-1.5)) + ((mod(t,4) >= 2 & mod(t,4) < 3) * (1.5));
tau = @(t) (0);


% Initial conditions (start from rest)
theta0 = pi/2;             % Initial angle (rad)
theta_dot0 = 0;         % Initial angular velocity (rad/s)
initial_conditions = [theta0; theta_dot0];

% Time span for the simulation
t_span = [0, 10];       % Simulate from t=0 to t=20 seconds (slower simulation)

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
    
    % Applied torque (from the smoothing function)
    torque = tau(t);  % Get the current torque value from the sigmoid function
    
    % Effective inertia
    I_eff = (2 / 3) * m_s * S^2 + m_t * S^2 * cos(theta1)^2;
    
    % Differential equation for theta1 (since theta2 = theta1 in the parallelogram)
    theta1_ddot = (spring_force_1 + spring_force_2 - torque) / I_eff;
    
    % Return the derivatives
    dydt = [theta_dot1; theta1_ddot];
end

% Solve the ODE using ode45
[t, y] = ode45(@(t, y) system_eqns(t, y, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau), t_span, initial_conditions);


% Define the function to calculate net torque (torque - spring forces)
function net_torque = equilibrium_func(theta, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau)
    % Calculate spring lengths at given theta
    l1 = sqrt((a * sin(theta))^2 + (a * cos(theta) + Q1)^2); % Spring 1 length
    l2 = sqrt((b * sin(theta))^2 + ((Q2 - T) + b * cos(theta))^2); % Spring 2 length
    
    % Derivatives of spring lengths
    dl1_dtheta = (a * cos(theta) + Q1) / l1;
    dl2_dtheta = ((Q2 - T) + b * cos(theta)) / l2;
    
    % Spring forces (passive forces, only reacting to displacement)
    spring_force_1 = k1 * (l1 - l1_o) * dl1_dtheta;
    spring_force_2 = k2 * (l2 - l2_o) * dl2_dtheta;
    
    % Applied torque (from the smoothing function)
    torque = tau(theta);  % Get the current torque value from the sigmoid function
    
    % Net torque (difference between applied torque and spring forces)
    net_torque = torque - (spring_force_1 + spring_force_2);
end

% Use fzero to find equilibrium points numerically (root of net_torque function)
theta_eq1 = fzero(@(theta) equilibrium_func(theta, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau), 0);  % Starting guess at 0
theta_eq2 = fzero(@(theta) equilibrium_func(theta, S, m_s, m_t, k1, k2, l1_o, l2_o, a, b, Q1, Q2, T, D, tau), pi/2);  % Starting guess at 90 degrees

disp(['Equilibrium at theta = 0: ', num2str(theta_eq1)]);
disp(['Equilibrium at theta = 90: ', num2str(theta_eq2)]);




theta_values = y(:,1); % in radians

% Pre-compute spring forces for each time instant
numPoints = length(t);
spring_force1 = zeros(numPoints,1);
spring_force2 = zeros(numPoints,1);
for i = 1:numPoints
    theta1 = theta_values(i);
    l1 = sqrt((a * sin(theta1))^2 + (a * cos(theta1) + Q1)^2);
    dl1_dtheta = (a * cos(theta1) + Q1) / l1;
    spring_force1(i) = k1 * (l1 - l1_o) * dl1_dtheta;
    
    l2 = sqrt((b * sin(theta1))^2 + ((Q2 - T) + b * cos(theta1))^2);
    dl2_dtheta = ((Q2 - T) + b * cos(theta1)) / l2;
    spring_force2(i) = k2 * (l2 - l2_o) * dl2_dtheta;
end

% Applied torque values over time (computed from the square wave)
torque_values = arrayfun(tau, t);

% Create a figure to plot the time-series data
figure;
subplot(4,1,1);
plot(t, theta_values, 'b-', 'LineWidth', 2);
ylabel('\theta (rad)');
title('Angle vs. Time');
grid on;

subplot(4,1,2);
plot(t, torque_values, 'm-', 'LineWidth', 2);
ylabel('Torque (N.m)');
title('Applied Torque vs. Time');
grid on;

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

% Plot the torque indication as text
h_data = text(0.1, 0.18, '', 'FontSize', 12, 'Color', 'm');  % Placeholder text for torque

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
    
    % Update displayed data (convert theta from rad to deg)
    currentThetaDeg = rad2deg(theta1);
    currentTorque = tau(t(i));
    currentSpringForce1 = spring_force1(i);
    currentSpringForce2 = spring_force2(i);
    
    data_str = sprintf('Theta: %.1f deg\nTorque: %.2f N.m\nSpring1: %.2f N\nSpring2: %.2f N', ...
        currentThetaDeg, currentTorque, currentSpringForce1, currentSpringForce2);
    set(h_data, 'String', data_str);
    
    % Pause for animation effect (slower simulation)
    pause(0.1);  % Slower animation, adjust this for desired speed
end

