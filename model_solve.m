%% Parameter Sweep to Find Bistable Equilibria
% Clear workspace and command window
clear;
clc;

%% Fixed Mechanism & Material Parameters
PLA_den = 1240;             % PLA density, kg/m^3
S = 0.05;                   % Side link length (m)
T = 0.085;                  % Top link length (m)
D = T;                      % Distance between pivots (m)
w_s = 0.005;                % Width of side link (m)
w_t = 0.0108;               % Width of top link (m)
m_s = PLA_den * w_s * S;    % Mass of each side link (kg)
m_t = PLA_den * w_t * T;    % Mass of the top link (kg)

% Fixed attachment (pivot) parameters
% These remain unchanged during the sweep:
a = 0.025;                % Distance from left pivot to spring 1 attachment (m)
b = 0.025;                % Distance from right pivot to spring 2 attachment (m)
Q1 = -0.03;               % Spring 1 fixed (bottom) attachment (m)
Q2 = D + 0.015;           % Spring 2 fixed (bottom) attachment (m)
E = 55158;                % Elastic modulus (N/m)

%% Define Sweep Ranges for Spring Geometry & Rest Lengths
% These are the design parameters that govern the computed spring stiffnesses.
% (Here we use a very coarse grid with two values per parameter. Increase the number
% of points if you wish to explore more combinations.)

p1_range = [0.1, 0.5];       % L spring1 (m)
q1_range = [0.05, 0.2];       % W spring1 (m)
r1_range = [0.01, 0.05];      % H spring1 (m)

p2_range = [0.1, 0.5];        % L spring2 (m)
q2_range = [0.05, 0.2];       % W spring2 (m)
r2_range = [0.01, 0.05];      % H spring2 (m)

l1_o_range = [0.1, 0.5];      % Rest length of spring 1 (m)
l2_o_range = [0.1, 0.5];      % Rest length of spring 2 (m)

%% Set Targets and Tolerances for Equilibrium
target_theta_0 = 0;         % target for tau_const = 0 [rad]
target_theta_90 = pi/2;     % target for tau_const = -1 [rad]
tol0 = 0.1;   % tolerance in radians around target 0 equilibrium
tol90 = 0.2;  % tolerance in radians around target 90 deg equilibrium

%% Prepare to store results
% We will store: p1, q1, r1, p2, q2, r2, l1_o, l2_o, theta_eq1, theta_eq2, err0, err90.
results = [];
count = 0;

%% Parameter Sweep Loops (Nested Loops)
for p1 = p1_range
    for q1 = q1_range
        for r1 = r1_range
            for p2 = p2_range
                for q2 = q2_range
                    for r2 = r2_range
                        for l1_o = l1_o_range
                            for l2_o = l2_o_range
                                % Compute spring stiffnesses
                                k1 = E * q1 * r1 / p1;  % Spring 1 stiffness (N/m)
                                k2 = E * q2 * r2 / p2;  % Spring 2 stiffness (N/m)
                                
                                % Define modular functions for spring lengths and their derivatives
                                l1_fun = @(theta) sqrt((a * sin(theta))^2 + (a * cos(theta) + Q1)^2);
                                dl1_fun = @(theta) (a * cos(theta) + Q1) / l1_fun(theta);
                                
                                l2_fun = @(theta) sqrt((b * sin(theta))^2 + ((Q2 - T) + b * cos(theta))^2);
                                dl2_fun = @(theta) ((Q2 - T) + b * cos(theta)) / l2_fun(theta);
                                
                                % Define constant torques for equilibrium search:
                                tau_const1 = 0;    % For equilibrium near theta = 0
                                tau_const2 = -1;   % For equilibrium near theta = 90 deg
                                
                                % Use fzero to find equilibrium angles:
                                % Initial guesses: 0 for the zero-torque condition,
                                % and pi/2 for the -1 torque condition.
                                try
                                    theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, ...
                                        l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const1), 0);
                                catch
                                    theta_eq1 = NaN;
                                end
                                try
                                    theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, ...
                                        l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const2), pi/2);
                                catch
                                    theta_eq2 = NaN;
                                end
                                
                                % Calculate errors relative to targets
                                err0 = abs(theta_eq1 - target_theta_0);
                                err90 = abs(theta_eq2 - target_theta_90);
                                
                                count = count + 1;
                                results(count,:) = [p1, q1, r1, p2, q2, r2, l1_o, l2_o, theta_eq1, theta_eq2, err0, err90];
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Convert Results to Table and Display
% Define variable names for clarity
varNames = {'p1','q1','r1','p2','q2','r2','l1_o','l2_o','theta_eq_tau0','theta_eq_tauNeg1','err0','err90'};
resultsTable = array2table(results, 'VariableNames', varNames);

% Display full table (or you can filter for the best candidates)
disp('Full Sweep Results:');
disp(resultsTable);

%% Optionally, filter results that meet the target tolerances
sel = results(results(:,11) < tol0 & results(:,12) < tol90, :);
if ~isempty(sel)
    disp('Parameter combinations meeting the target equilibria conditions (within tolerance):');
    selTable = array2table(sel, 'VariableNames', varNames);
    disp(selTable);
else
    disp('No parameter combination met the target conditions within the specified tolerances.');
end

%% --- Local Function Definitions ---
function netTorque = equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const)
    % equilibrium_func computes the net torque (applied constant torque minus sum of spring torques)
    % for the given configuration angle theta.
    % Equilibrium is when netTorque == 0.
    l1_val = l1_fun(theta);
    l2_val = l2_fun(theta);
    dl1_val = dl1_fun(theta);
    dl2_val = dl2_fun(theta);
    
    springForce1 = k1 * (l1_val - l1_o) * dl1_val;
    springForce2 = k2 * (l2_val - l2_o) * dl2_val;
    
    netTorque = tau_const - (springForce1 + springForce2);
end
