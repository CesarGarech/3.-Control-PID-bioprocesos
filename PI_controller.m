function F = PI_controller(t, S, S_setpoint, params)
    persistent integral
    if isempty(integral)
        integral = 0;
    end
    
    % PI parameters
    Kp = params.Kp;
    Ki = params.Ki;
    
    % Error
    error = S_setpoint - S;
    
    % Integral term
    integral = integral + error * params.Ts;
    
    % PI control law
    F = Kp * error + Ki * integral;
    
    % Ensure the feed rate is non-negative
    F = max(F, 0);
end
