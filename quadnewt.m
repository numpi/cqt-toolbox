function X = quadnewt(Am1, A0, A1)
%QUADNEWT Solve A1 X^2 + A0 X + Am1 = 0

debug = false;

maxit = 10;

tol = sqrt( cqtoption('threshold') );

% Compute the residual -- we save S because it will be useful later on, and
% in this way we can perform the first iteration right away. 
X = Am1;

for j = 1 : maxit
    
    % Compute the residual -- we save S because it will be useful later on
    S = (A1 * X + A0); iS = inv(S);
    R = S * X + Am1;
    
    nrmR = norm(R, inf);
    
    if debug
        fprintf('Newton :: it %d, res = %e\n', j, nrmR);
    end
    
    if nrmR < tol
        break;
    end
    
    % Solve the equation for the Newton iteration   
    M = iS * A1;
    gamma = sqrt( norm(M, inf) / norm(X, inf) );
    
    % Here we use the fact that iS * R = X + iS * Am1;
    D = cqtstein(M / gamma, X * gamma, X + iS * Am1);
    
    X = X + D;
end

end

