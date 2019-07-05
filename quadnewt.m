function X = quadnewt(Am1, A0, A1, varargin)
%QUADNEWT Solve A1 X^2 + A0 X + Am1 = 0
%
% X = QUADNEWT(AM1, A0, A1) solves the quadratic matrix equation
%
%              AM1 + A0 * X + A1 * X^2 = 0            (1), 
%
%     using the Newton iteration. It is assumed that, at every step of the
%     Newton method, the Stein equatino MXN + X + R = 0 arising can be
%     solved assuming that the spectral radii of M and N are smaller than
%     1. This is guaranteed under certain assumption if the coefficients
%     arise from a 2-level QBD process. 
%
% X = QUADNEWT(AM1, A0, A1, KEY1, VAL1, ..., KEY_N, VAL_N) allows to
%     specify additional options in key / value pair. The available options
%     are:
%      - 'tol': Set the tolerance to stop the Newton iteration. The default
%               is 1e2 * cqtoption('threshold'); 
%      - 'debug': Can be true or false; if true, shows the residuals during
%                 the iterations. The default is false. 
%      - 'maxit': Maximum number of allowed Newton iterations. The default
%                 is 15. 
%      - 'X': An optional input parameter -- if not given, the iteration
%             starts from the zero matrix. 

% parse Inputs
p = inputParser;

addParameter(p, 'debug', false);
addParameter(p, 'tol', 1e1 * ( cqtoption('threshold') ));
addParameter(p, 'maxit', 100);
addParameter(p, 'X', -A0 \ Am1);
addParameter(p, 'method', 'galerkin');

parse(p, varargin{:});

debug = p.Results.debug;
maxit = p.Results.maxit;
tol   = p.Results.tol;
method = p.Results.method;

X = p.Results.X;

j = 1;

while j < maxit
    
    % Compute the residual -- we save S because it will be useful later on
    S = (A1 * X + A0); iS = inv(S);
    R = S * X + Am1;
    
    nrmR = norm(R, inf);
    
    [rm, rp] = symbol(R);
    
    % If the symbol
    if norm(rm, 1) + norm(rp, 1) < cqtoption('threshold')
        [U, V] = correction(R);
        R = cqt([], [], U, V);
    end
    
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
    D = cqtstein(M / gamma, X * gamma, X + iS * Am1, 'tol', ...
        min(sqrt(tol), max(tol / 10, nrmR^2)), 'method', method, 'debug', debug);
    
    X = X + D;
    
    j = j + 1;
end

end

