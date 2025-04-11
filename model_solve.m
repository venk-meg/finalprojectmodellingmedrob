%% Parameters
PLA_den = 1240;         % PLA density kg/m3
S = 0.05;               % Length of each side link (m)
T = 0.085;              % Length of the top link (m)
w_s = 0.005;            % Width of side link (m)
w_t = 0.0102;           % Width of top link (m)
m_s = PLA_den * w_s * S;   % Mass of each side link (kg)
m_t = PLA_den * w_t * T;   % Mass of the top link (kg)
k1 = 100;               % Spring 1 stiffness (N/m)
k2 = 0;                 % Spring 2 stiffness (N/m) (if unused, set to zero)
l1_o = 0.03;            % Rest length of spring 1 (m)
l2_o = 0.03;            % Rest length of spring 2 (m)
a = 0.025;              % Distance along OA from O to the attachment point "a"
b = 0.025;              % Distance along CB from C to the attachment point "b"
D = 0.085;              % Distance between fixed pivots O and C (m)
Q1 = -0.05;             % Location of spring 1’s bottom attachment (m) (x-coordinate)
Q2 = D + 0.09;          % Location of spring 2’s bottom attachment (m) (x-coordinate)

%% Define the geometry functions for springs

% l1: length of spring 1 as a function of theta
l1_fun = @(theta) sqrt((a * sin(theta))^2 + (a * cos(theta) + Q1)^2);
% Derivative dl1/dtheta:
dl1_dtheta_fun = @(theta) (a * cos(theta) + Q1) / l1_fun(theta);

% l2: length of spring 2 as a function of theta
l2_fun = @(theta) sqrt((b * sin(theta))^2 + ((Q2 - T) + b * cos(theta))^2);
% Derivative dl2/dtheta:
dl2_dtheta_fun = @(theta) ((Q2 - T) + b * cos(theta)) / l2_fun(theta);

%% Define the equilibrium function
% Equilibrium (net torque) function:
% f(theta) = spring force contributions - applied torque (tau_val)
% Here, applied torque is taken as a constant tau_val.
eqfun = @(theta, tau_val) k1*(l1_fun(theta) - l1_o)*dl1_dtheta_fun(theta) + ...
                           k2*(l2_fun(theta) - l2_o)*dl2_dtheta_fun(theta) - tau_val;
                       
%% Find Equilibrium Points Using fzero
% For equilibrium when tau_val = 0 (expected equilibrium at theta = 0)
tau_val1 = 0;
theta_guess1 = 0;  % initial guess near 0
theta_eq1 = fzero(@(theta) eqfun(theta, tau_val1), theta_guess1);

% For equilibrium when tau_val = 1 (you may also search for tau_val = -1)
tau_val2 = 1;
theta_guess2 = pi/2;  % initial guess near 90 degrees (pi/2 rad)
theta_eq2 = fzero(@(theta) eqfun(theta, tau_val2), theta_guess2);

fprintf('Equilibrium with tau = %.2f N.m: theta = %.4f rad (%.2f deg)\n', tau_val1, theta_eq1, rad2deg(theta_eq1));
fprintf('Equilibrium with tau = %.2f N.m: theta = %.4f rad (%.2f deg)\n', tau_val2, theta_eq2, rad2deg(theta_eq2));

%% Optionally, plot the equilibrium function vs. theta to visualize the roots
theta_vals = linspace(0, pi/2, 200);
f1_vals = arrayfun(@(theta) eqfun(theta, tau_val1), theta_vals);
f2_vals = arrayfun(@(theta) eqfun(theta, tau_val2), theta_vals);

figure;
subplot(2,1,1);
plot(theta_vals, f1_vals, 'b-', 'LineWidth', 2);
hold on;
plot(theta_eq1, eqfun(theta_eq1, tau_val1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('\theta (rad)'); ylabel('Net torque (N.m)');
title(['Equilibrium for \tau = ' num2str(tau_val1) ' N.m']);
grid on;

subplot(2,1,2);
plot(theta_vals, f2_vals, 'b-', 'LineWidth', 2);
hold on;
plot(theta_eq2, eqfun(theta_eq2, tau_val2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('\theta (rad)'); ylabel('Net torque (N.m)');
title(['Equilibrium for \tau = ' num2str(tau_val2) ' N.m']);
grid on;
