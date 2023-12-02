clear
close all
clc

%% Initial Conditions
X0 = 0.01;             
S0 = 20.0;             
P0 = 0.0;             
V0 = 1.0;             
initial_conditions = [X0; S0; P0; V0];

%% Control opt.
S_setpoint = 18.0;     
tspan = [0 40];       
params.Kp = 10.5;
params.Ki = 0.001;
params.Ts = 0.05;     
options = odeset('NonNegative', 1:4);

%% EKF Initialization
mu_max_est = 0.4;  
Y_XS_est = 0.6;  
alpha_est = 0.1;   
beta_est = 0.05;   
params_est = [mu_max_est; Y_XS_est; alpha_est; beta_est];

X_est = [X0; S0; P0; V0; params_est];
P0_est = eye(8)*0.1;
Q = eye(8)*0.001;
R = 0.01;

%% Simulation loop with EKF
time_points = tspan(1):params.Ts:tspan(2);
num_points = numel(time_points);
S_values = zeros(1, num_points);
X_values = zeros(1, num_points);
F_values = zeros(1, num_points);
X_est_values = zeros(length(X_est), num_points);

for i = 1:num_points
    t = time_points(i);
    
    if i == 1
        Y_current = initial_conditions;
    else
        [~, Y] = ode45(@(t,Y) bioreactor_model(t, Y, F_values(i-1), X_est(5:8)), [time_points(i-1), t], Y_current, options);
        Y_current = Y(end, :);
    end
    
    S_current = Y_current(2);
    X_current = Y_current(1);
    S_values(i) = S_current;
    X_values(i) = X_current;
    
    if i > 1
        z = X_current + sqrt(R)*randn;
        [X_est, P0_est] = ekf_predict(X_est, P0_est, Q, @process_model);
        [X_est, P0_est] = ekf_update(X_est, P0_est, z, R, @measurement_model);
    end
    
    X_est_values(:, i) = X_est;
    
    if Y_current(4) < 10
        F_values(i) = PI_controller(t, S_current, S_setpoint, params);
    else
        F_values(i) = 0;
    end
end

%% Plot the results
% ... (Your existing plotting code)
% ...

%% Functions for EKF
function X_next = process_model(X)
    % Simple example; replace with your bioreactor equations
    mu = X(5) * X(2) / (0.1 + X(2));
    X_next = X + [mu * X(1); -mu * X(1)/X(6); X(7) * mu * X(1); 0; 0; 0; 0; 0];
end

function Hx = measurement_model(X)
    Hx = X(1);
end

function [X_pred, P_pred] = ekf_predict(X_est, P_est, Q, process_model)
    X_pred = process_model(X_est);
    A = jacobian(@(X) process_model(X), X_est);
    P_pred = A * P_est * A' + Q;
end

function [X_upd, P_upd] = ekf_update(X_est, P_est, z, R, measurement_model)
    Hx = measurement_model(X_est);
    H = jacobian(@(X) measurement_model(X), X_est);
    K = P_est * H' / (H * P_est * H' + R);
    X_upd = X_est + K * (z - Hx);
    P_upd = (eye(size(K,1)) - K * H) * P_est;
end

function J = jacobian(fun, X)
    n = length(X);
    f = fun(X);
    J = zeros(length(f), n);
    h = n * eps;
    for i = 1:n
        X_temp = X;
        X_temp(i) = X_temp(i) + h;
        f_temp = fun(X_temp);
        J(:, i) = (f_temp - f) / h;
    end
end

function dY = bioreactor_model(t, Y, F, params)
    % Your bioreactor model goes here
    % Example:
    mu = params(1) * Y(2) / (0.1 + Y(2));
    dY = [mu * Y(1); -mu * Y(1)/params(2); params(3) * mu * Y(1); F];
end

function F = PI_controller(t, S_current, S_setpoint, params)
    % Your existing PI_controller function
    % Example:
    e = S_setpoint - S_current;
    F = params.Kp * e + params.Ki * e * params.Ts;
end
