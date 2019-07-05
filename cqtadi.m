function X = cqtadi(A, B, C, alpha, beta, varargin)
%CQTADI CQT-version of the ADI iteration

p = inputParser;

addParameter(p, 'tol', cqtoption('threshold'));
addParameter(p, 'maxit', 1000);
addParameter(p, 'debug', false);

parse(p, varargin{:});

tol   = max(p.Results.tol, 10 * cqtoption('threshold'));
maxit = p.Results.maxit;
debug = p.Results.debug;

it = 1;

c = 1;

X = cqt([], []);
I = cqt(1, 1);

res = inf;

while it < maxit
    Xold = X;
    X = - (A - beta(c) * I) \ ( X * (B + beta(c) * I) + C );
    X = - ((A - alpha(c) * I) * X + C) / (B + alpha(c) * I);
    
    oldres = res;
    res = norm(X - Xold, inf); % / norm(X, inf);
    
    if debug
        fprintf('ADI :: Iteration %d, res = %e, %e\n', it, res, norm(A*X + X*B + C, inf));
    end
    
    if res < tol % || oldres < res
        break;
    end
    
    c = 1 + mod(c, length(alpha));
    it = it + 1;
end

end

