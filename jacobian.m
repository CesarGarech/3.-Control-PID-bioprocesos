function J = jacobian(fun, X)
% NUMERICAL_JACOBIAN computes the Jacobian matrix of the function fun at the point X
    % using central difference method.
    %
    % Inputs:
    %    fun - Function handle
    %    X   - Point at which the Jacobian is to be computed
    %
    % Output:
    %    J   - Jacobian matrix
    
    n = length(X);
    f = fun(X);
    J = zeros(length(f), n);
    h = sqrt(eps) * (1 + abs(X));
    
    for i = 1:n
        X_temp1 = X;
        X_temp2 = X;
        X_temp1(i) = X_temp1(i) + h(i);
        X_temp2(i) = X_temp2(i) - h(i);
        f_temp1 = fun(X_temp1);
        f_temp2 = fun(X_temp2);
        J(:, i) = (f_temp1 - f_temp2) / (2 * h(i));
    end
end