clear all
close all
clc

%% Initial Conditions for the model
X0 = 0.01;             % g/L
S0 = 20.0;             % g/L
P0 = 0.0;              % g/L
V0 = 1.0;              % L
initial_conditions = [X0; S0; P0; V0];

%% Control options.
% Setpoint
S_setpoint = 18.0;     % g/L

% PI controller parameters
params.Kp = 10.5;
params.Ki = 0.001;
params.Ts = 0.1;      % Sampling time, h

%% EKF Initialization
mu_max_est = 0.83;
Y_XS_est = 0.8;
alpha_est = 0.05;
beta_est = 0.002;
params_est = [mu_max_est; Y_XS_est; alpha_est; beta_est];

initialState = [X0; S0; P0; V0; params_est];
P0_est = eye(8).*[1e-3 1e3 1e-3 1e-1 1e-1 1e-1 1e-1 1e-1];
Q = eye(8).*[1e2 1e0 1e-3 1e-3 1e-1 1e-1 1e-1 1e-1];
R = 0.0;

% Create the trackingEKF object
ekf = trackingEKF(@process_modelDT, @measurement_model, ...
    'State', initialState, 'StateCovariance', P0_est, ...
    'ProcessNoise', Q, 'MeasurementNoise', R);

%% Simulation loop
% Simulation Time
tspan = [0 40];       % h
time_points = tspan(1):params.Ts:tspan(2);
num_points = numel(time_points);
S_values = zeros(1, num_points);
X_values = zeros(1, num_points);
P_values = zeros(1, num_points);
V_values = zeros(1, num_points);
F_values = zeros(1, num_points);
X_est_values = zeros(length(initialState), num_points);

for i = 1:num_points
    t = time_points(i);

    if i == 1
        Y_current = initial_conditions;
        F_current = 0; % Initial feed rate can be set to zero or an initial value
    else
        % ODE solver options
        options = odeset('NonNegative', 1:4);
        [~, Y] = ode15s(@(t,Y) bioreactor_model(t, Y, F_current, ekf.State(5:8)), [time_points(i-1), t], Y_current, options);
        Y_current = Y(end, :);
        % Get the previous feed rate
        F_current = F_values(i-1);
    end

    % Store the actual states
    S_values(i) = Y_current(2);
    X_values(i) = Y_current(1);
    P_values(i) = Y_current(3);
    V_values(i) = Y_current(4);

    % Simulate the measurement
    z = X_values(i) + sqrt(R)*randn;

    % EKF Predict and Update
    predict(ekf, F_current);
    correct(ekf, z);

    % Store the estimated states
    X_est = ekf.State;
    X_est_values(:, i) = X_est;

    S_estim = X_est(2);
    if V_values(i) < 10
        F_values(i) = PI_controller(t, S_estim, S_setpoint, params);
    else
        F_values(i) = 0;
    end
    F_current = F_values(i);
end

SetP = ones(length(S_values),1) * S_setpoint;

%% Plot the results
figure(1);
subplot(2,2,1)
plot(time_points, S_values, 'LineWidth', 2);hold on;plot(time_points,SetP,'LineWidth', 2);
xlabel('Time (h)')
ylabel('Substrate Concentration (g/L)')
title('Substrate Concentration vs Time')

subplot(2,2,2)
plot(time_points, X_values, 'r', 'LineWidth', 2);
xlabel('Time (h)')
ylabel('Biomass Concentration')
title('Biomass concentration vs Time')

subplot(2,2,3)
plot(time_points, P_values, 'r', 'LineWidth', 2);
xlabel('Time (h)')
ylabel('Product Concentration')
title('Product concentration vs Time')

sgtitle('Fed-Batch Bioreactor Control')

subplot(2,2,4)
plot(time_points, F_values, 'r', 'LineWidth', 2);
xlabel('Time (h)')
ylabel('Feed Rate (L/h)')
title('Feed Rate vs Time')

figure(2);
subplot(2,2,1)
plot(time_points, X_values, 'LineWidth', 2);hold on;plot(time_points,X_est_values(1,:),'*','LineWidth', 2);
xlabel('Time (h)')
ylabel('Biomass Concentration (g/L)')
legend('Plant','EKF',Location='best')
title('Biomass Concentration vs Time')

subplot(2,2,2)
plot(time_points, S_values, 'LineWidth', 2);hold on;plot(time_points,X_est_values(2,:),'*','LineWidth', 2);
xlabel('Time (h)')
ylabel('Substrate Concentration (g/L)')
legend('Plant','EKF',Location='best')
title('Substrate Concentration vs Time')

subplot(2,2,3)
plot(time_points, P_values, 'LineWidth', 2);hold on;plot(time_points,X_est_values(3,:),'*','LineWidth', 2);
xlabel('Time (h)')
ylabel('Product Concentration (g/L)')
legend('Plant','EKF',Location='best')
title('Product Concentration vs Time')

subplot(2,2,4)
plot(time_points, V_values, 'LineWidth', 2);hold on;plot(time_points,X_est_values(4,:),'*','LineWidth', 2);
xlabel('Time (h)')
ylabel('Volume (L)')
legend('Plant','EKF',Location='best')
title('Volume vs Time')

sgtitle('EKF State-Estimation')

figure(3);
umax = ones(length(X_values),1)*initialState(5);
subplot(2,2,1)
plot(time_points,X_est_values(5,:),'*','LineWidth', 2);hold on;plot(time_points,umax,'LineWidth', 2);
xlabel('Time (h)')
ylabel('Umax (h^-1)')
legend('EKF','Plant',Location='best')
title('Umax vs Time')

yxs = ones(length(X_values),1)*initialState(6);
subplot(2,2,2)
plot(time_points,X_est_values(6,:),'*','LineWidth', 2);hold on;plot(time_points,yxs,'LineWidth', 2);
xlabel('Time (h)')
ylabel('Yxs (g/g)')
legend('EKF','Plant',Location='best')
title('Yxs vs Time')

alpha = ones(length(X_values),1)*initialState(7);
subplot(2,2,3)
plot(time_points,X_est_values(7,:),'*','LineWidth', 2);hold on;plot(time_points,alpha,'LineWidth', 2);
xlabel('Time (h)')
ylabel('alpha (Adim)')
legend('EKF','Plant',Location='best')
title('alpha vs Time')

beta = ones(length(X_values),1)*initialState(8);
subplot(2,2,4)
plot(time_points,X_est_values(8,:),'*','LineWidth', 2);hold on;plot(time_points,beta,'LineWidth', 2);
xlabel('Time (h)')
ylabel('beta (Adim)')
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%.4f';
legend('EKF','Plant',Location='best')
title('beta vs Time')

sgtitle('EKF Parameter-Estimation')