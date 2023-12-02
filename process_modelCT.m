function vard = process_modelCT(X, F)
    % Extract states and parameters from the input
    X_bio = X(1);
    S = X(2);
    P = X(3);
    V = X(4);
    mu_max = X(5);
    Y_XS = X(6);
    alpha = X(7);
    beta = X(8);
    
 
    % Monod Growth kinetics
    Ks = 0.08;          % g/L
    Smax =30; %g/L
    ms = 0.00001;
    mu = mu_max * S / (Ks + S) * (1 - (S / Smax));
    
    % External values, which could be additional inputs to the function
%     F = 20; % Flow rate, which comes from your controller
    S_in = 20; %g/L
    
    % Mass balances
    dXdt = mu * X_bio - (F / V) * X_bio;
    dSdt = (F / V) * (S_in - S) - ((mu / Y_XS) + ms) * X_bio;
    dPdt = (alpha * mu + beta) * X_bio - (F / V) * P;
    dVdt = F;
    
    % No dynamics for the parameters in the prediction step, they are considered constant
    dmu_max_dt = 0;
    dY_XS_dt = 0;
    dalpha_dt = 0;
    dbeta_dt = 0;

    vard = [dXdt; dSdt; dPdt; dVdt; dmu_max_dt; dY_XS_dt; dalpha_dt; dbeta_dt];
    
%     % Return the next state prediction
%     X_next = X + dt*[dXdt; dSdt; dPdt; dVdt; dmu_max_dt; dY_XS_dt; dalpha_dt; dbeta_dt];
end
