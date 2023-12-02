function [X_upd, P_upd] = ekf_update(X_est, P_est, z, R, measurement_model)
    Hx = measurement_model(X_est);
    H = jacobian(@(X) measurement_model(X), X_est);
    K = P_est * H' / (H * P_est * H' + R);
    X_upd = X_est + K * (z - Hx);
    P_upd = (eye(size(K,1)) - K * H) * P_est;
end