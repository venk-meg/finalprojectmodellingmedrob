%% Parameter Sweep for Spring and Attachment Parameters with Filtering
% Clear workspace and command window
clear;
clc;

%% Fixed System Parameters
E = 55158;       % Elastic modulus (N/m)
T = 0.085;       % Length of the top link (m)
D = T;           % Distance between pivots (m)

%% Define Ranges for Spring Parameters
% For demonstration, we use two values for each parameter (expand as needed)

Q_1_vals = [-0.03];     
Q_2_vals = [D + 0.03]; 

p1_vals = [0.01,0.03, 0.05];      % L spring1 (m)
q1_vals = [0.01, 0.02];     % W spring1 (m)
r1_vals = [0.001, 0.003, 0.005];    % H spring1 (m)

p2_vals = p1_vals;     % L spring2 (m)
q2_vals = q1_vals;     % W spring2 (m)
r2_vals = r1_vals;      % H spring2 (m)

l1_o_vals = [0.01, 0.03, 0.05];    % Rest length of spring 1 (m)
l2_o_vals = l1_o_vals;    % Rest length of spring 2 (m)

a_vals = [0.01, 0.03, 0.05];     % Distance from left pivot to spring 1 attachment (m)
b_vals = a_vals;     % Distance from right pivot to spring 2 attachment (m)

%% Initialize Results Storage
% Columns: [p1, q1, r1, p2, q2, r2, l1_o, l2_o, a, b, ...
%           theta_eq1_rad, theta_eq1_deg, theta_eq2_rad, theta_eq2_deg]
results = [];

%% Nested Loops Over Parameter Combinations
for Q_1 = Q_1_vals
    for Q_2 = Q_2_vals
        for p1 = p1_vals
            for q1 = q1_vals
                for r1 = r1_vals
                    for p2 = p2_vals
                        for q2 = q2_vals
                            for r2 = r2_vals
                                for l1_o = l1_o_vals
                                    for l2_o = l2_o_vals
                                        for a = a_vals
                                            for b = b_vals
                                                % Compute spring stiffnesses for the current combination
                                                k1 = E * q1 * r1 / p1;
                                                k2 = E * q2 * r2 / p2;
                                                
                                                % Define modular spring functions for this combination
                                                l1_fun = @(theta) sqrt((a * sin(theta))^2 + (a * cos(theta) + Q_1)^2);
                                                dl1_fun = @(theta) (a * cos(theta) + Q_1) / l1_fun(theta);
                                                
                                                l2_fun = @(theta) sqrt((b * sin(theta))^2 + ((Q_2 - T) + b * cos(theta))^2);
                                                dl2_fun = @(theta) ((Q_2 - T) + b * cos(theta)) / l2_fun(theta);
                                                
                                                % Find equilibrium for tau_const = 0 (target near θ = 0)
                                                try
                                                    theta_eq1 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, ...
                                                        l1_fun, dl1_fun, l2_fun, dl2_fun, 0), 0);
                                                catch
                                                    theta_eq1 = NaN;
                                                end
                                                
                                                % Find equilibrium for tau_const = -1 (target near θ = 90° or pi/2 rad)
                                                try
                                                    theta_eq2 = fzero(@(theta) equilibrium_func(theta, k1, k2, l1_o, l2_o, ...
                                                        l1_fun, dl1_fun, l2_fun, dl2_fun, -1), pi/2);
                                                catch
                                                    theta_eq2 = NaN;
                                                end
                                                
                                                % Append the current parameter set and computed equilibrium angles
                                                results = [results; p1, q1, r1, p2, q2, r2, l1_o, l2_o, a, b, Q_1, Q_2,...
                                                    theta_eq1, rad2deg(theta_eq1), theta_eq2, rad2deg(theta_eq2)];
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
    end
end

%% Create a Table for All Results
varNames = {'p1','q1','r1','p2','q2','r2','l1_o','l2_o','a','b', 'Q_1', 'Q_2', ...
    'theta_eq1_rad','theta_eq1_deg','theta_eq2_rad','theta_eq2_deg'};
resultsTable = array2table(results, 'VariableNames', varNames);

%% Apply Filtering Based on Thresholds
% Define thresholds in degrees:
% For tau_const = 0, we want equilibrium near 0 deg (absolute value less than threshold0)
% For tau_const = -1, we want equilibrium near 90 deg (within threshold90 of 90)
threshold0 = 7;   % degrees threshold for equilibrium near 0
threshold90 = 7;  % degrees threshold for equilibrium near 90

% Logical indices for rows that satisfy both conditions:
filterIdx = (abs(results(:,12)) < threshold0) & (abs(results(:,14) - 90) < threshold90);
filteredResults = results(filterIdx, :);
filteredTable = array2table(filteredResults, 'VariableNames', varNames);

%% Display the Full and Filtered Results
disp('All Results:');
disp(resultsTable);

disp('Filtered Results (Equilibrium near 0 within ±10° and near 90° within ±10°):');
disp(filteredTable);

%% Local Function Definition
function netTorque = equilibrium_func(theta, k1, k2, l1_o, l2_o, l1_fun, dl1_fun, l2_fun, dl2_fun, tau_const)
    % Compute spring lengths and their derivatives
    l1 = l1_fun(theta);
    l2 = l2_fun(theta);
    dl1 = dl1_fun(theta);
    dl2 = dl2_fun(theta);
    
    % Compute spring forces: F = k*(l - l0) * (dl/dtheta)
    springForce1 = k1 * (l1 - l1_o) * dl1;
    springForce2 = k2 * (l2 - l2_o) * dl2;
    
    % Equilibrium condition: tau_const is balanced by the spring forces
    netTorque = tau_const - (springForce1 + springForce2);
end
