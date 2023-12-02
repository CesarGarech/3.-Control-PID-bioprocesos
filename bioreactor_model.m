function dYdt = bioreactor_model(t, Y, F, Params)
    
    % Parameters
    mu_max = Params(1);      % 1/h
    Y_XS = Params(2);        % g/g
    Ks = 0.08;          % g/L
    S_in = 20;         % g/L (substrate concentration in feed)
    Smax = 30;         %g/L
    ms = 0.00001;
    alfa = Params(3);
    beta = Params(4);

    % Current states
    X = Y(1);          % Biomass concentration (g/L)
    S = Y(2);          % Substrate concentration (g/L)
    P = Y(3);          % Product concentration (g/L)
    V = Y(4);          % Volume (L)
    
    % Monod Growth kinetics
    mu = mu_max * S / (Ks + S)* (1-(S/(Smax)));
    
    % Mass balances
    dXdt = mu * X - (F/V) * X;
    dSdt = (F/V) * (S_in - S) - ((mu / Y_XS) + ms) * X;
    dPdt = (alfa * mu + beta) * X - (F/V) * P;
    dVdt = F;
    
    dYdt = [dXdt; dSdt; dPdt; dVdt];
end
