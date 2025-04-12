%% Parameter Sweep for Bistable Equilibria
% Clear workspace and command window
clear;
clc;

%% Global Geometry & Material Parameters (fixed)
PLA_den = 1240;        % PLA density kg/m3
S = 0.05;              % Length of each side link (m)
T = 0.085;             % Length of the top link (m)
D = T;                 % Distance between pivots (m)
w_s = 0.005;           % Width of side link (m)
w_t = 0.0108;          % Width of top link (m)
m_s = PLA_den * w_s * S;   % Mass of each side link (kg)
m_t = PLA_den * w_t * T;   % Mass of the top link (kg)

E = 55158;             % Elastic modulus (N/m)
% For simplicity, we fix these two (but you could sweep them too)
q1 = 0.2; r1 = 0.05;    % Dimensions for spring 1 (width and height)
q2 = 0.2; r2 = 0.05;    % Dimensions for spring 2 (width and height)

%% Define Ranges for the Spring Parameters to Sweep
% Here we choose only a subset for demonstration.
p1_range    = linspace(0.1, 0.5, 5);    % L spring1 (m)
p2_range    = p1_range;                 % Let spring2 be similar (can be swept independently)
l1o_range   = linspace(0.1, 0.5, 5);      % Rest length of spring 1 (m)
l2o_range   = l1o_range;                % Rest length of spring 2 (m)
a_range     = linspace(0.02, 0.05, 5);    % Attachment point on left side (m)
b_range     = a_range;                  % Attachment point on right side (m)
Q1_range    = linspace(0, -0.1, 5);       % Fixed bottom attachment for spring 1 (from 0 to -0.1 m)
Q2_range    = linspace(D, D+0.1, 5);      % Fixed bottom attachment for spring 2 (from D to D+0.1 m)

%% Equilibrium Target Conditions and Tolerances
% We want:
%   For τ_const = 0, equilibrium near θ = 0 rad.
%   For τ_const = -1, equilibrium near θ = π/2 rad.
tol_eq_low = 0.1;    % tolerance for equilibrium near 0 (in rad)
tol_eq_high = 0.1;   % tolerance for equilibrium near π/2 (in rad)

%% Preallocate a results array to save matching parameter sets
% We'll store: [p1, l1_o, a, Q1, p2, l2_o, b, Q2, theta_eq1, theta_eq2]
results = [];

%% Parameter Sweep Loop
fprintf('Starting parameter sweep...\n');
for p1 = p1_range
    for p2 = p2_range
        % Compute stiffness for each spring (they will be recalculated once we choose p1 and p2)
        k1 = E * q1 * r1 / p1;
        k2 = E * q2 * r2 / p2;
        for l1_o = l1o_range
            for l2_o = l2o_range
                for a = a_range
                    for b = b_range
                        for Q1 = Q1_range
                            for Q2 = Q2_range
                                % Define the modular functions for spring 1 and spring 2 given current parameters:
                                l1_fun = @(theta) sqrt((a * sin(theta))^2 + (a * cos(theta) + Q1)^2);
                                dl1_fun = @(theta) (a * cos(theta) + Q1) / l1_fun(theta);
                                
                                l2_fun = @(theta) sqrt((b * sin(theta))^2 + ((Q2 - T) + b * cos(theta))^2);
                                dl2_fun = @(theta) ((Q2 - T) + b * cos(theta)) / l2_fun(theta);
                                
                                % Find equilibrium for tau_const = 0 (target: theta near 0)
                                try
                                    theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, 0), 0);
                                catch
                                    theta_eq1 = NaN;
                                end
                                % Find equilibrium for tau_const = -1 (target: theta near pi/2)
                                try
                                    theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, -1), pi/2);
                                catch
                                    theta_eq2 = NaN;
                                end
                                
                                % Check if both solutions are real numbers and meet the targets.
                                % For tau_const = 0, we want |theta_eq1| < tol_eq_low.
                                % For tau_const = -1, we want |theta_eq2 - (pi/2)| < tol_eq_high.
                                if ~isnan(theta_eq1) && ~isnan(theta_eq2)
                                    if abs(theta_eq1) < tol_eq_low && abs(theta_eq2 - (pi/2)) < tol_eq_high
                                        % Store the current parameter combination and the equilibrium values.
                                        results = [results; p1, l1_o, a, Q1, p2, l2_o, b, Q2, theta_eq1, theta_eq2];
                                        fprintf('Match found: p1=%.3f, l1_o=%.3f, a=%.3f, Q1=%.3f, p2=%.3f, l2_o=%.3f, b=%.3f, Q2=%.3f, theta_eq1=%.3f, theta_eq2=%.3f\n',...
                                            p1, l1_o, a, Q1, p2, l2_o, b, Q2, theta_eq1, theta_eq2);
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Display final results
if isempty(results)
    disp('No parameter combination met the desired criteria.');
else
    disp('The following parameter combinations produced the targeted equilibria:');
    % Columns: [p1, l1_o, a, Q1, p2, l2_o, b, Q2, theta_eq1, theta_eq2]
    disp(results);
end

%% Local Function Definitions
function netTorque = equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const)
    % Computes the net torque:
    % netTorque = tau_const - (springForce1 + springForce2)
    % An equilibrium is when netTorque = 0.
    l1 = l1_fun(theta);
    l2 = l2_fun(theta);
    dl1 = dl1_fun(theta);
    dl2 = dl2_fun(theta);
    
    springForce1 = k1 * (l1 - l1_o) * dl1;
    springForce2 = k2 * (l2 - l2_o) * dl2;
    
    netTorque = tau_const - (springForce1 + springForce2);
end
