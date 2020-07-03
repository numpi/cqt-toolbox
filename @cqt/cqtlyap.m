function X = cqtlyap(varargin)
%CQTLYAP Lyapunov and Sylvester Solver%
% 
% X = CQT(A, B, C) solves the Sylvester equation AX + XB + C = 0.

A = varargin{1};
B = varargin{2};
C = varargin{3};

p = inputParser;

addParameter(p, 'poles', []);
addParameter(p, 'tol', cqtoption('threshold'));
addParameter(p, 'debug', false);
addParameter(p, 'maxit', inf);

parse(p, varargin{4:end});

poles = p.Results.poles;
tol   = p.Results.tol;
debug = p.Results.debug;
maxit = p.Results.maxit;

if ~isinf(max(size(A)))
    error('cqtlyap is only supported for infinite matrices');
end

% Solve the equation on the Toeplitz part by evaluation / interpolation
[am, ap] = symbol(A);
[bm, bp] = symbol(B);
[cm, cp] = symbol(C);

% Construct the symbol of X
if isempty(cm) && isempty(cp)
    xm = [];
    xp = [];
else
    [xm, xp] = evinterp(@(a,b,c) c ./ (a + b), tol, am, ap, bm, bp, -cm, -cp);
end

% Solve the correction equation
X = cqt(xm, xp, [], [], size(A, 2), size(B, 1));
R = A*X + X*B + C;

[RU, RV] = correction(R);

% Solve the equation for the correction
if isempty(poles)
    [Xu, Xv] = ek_sylv(A, B, RU, RV, 2 * maxit * size(RU, 2), ...
        @(r, n) r < tol, debug, 2);
else
    [Xu, Xv] = rk_sylv(poles, A, B, RU, RV, maxit * size(RU, 2), ...
        @(r, n) r < tol, debug, 2);
end

X = X + cqt([], [], Xu, Xv);

end

