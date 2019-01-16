function X = cqtlyap(varargin)
%CQTLYAP Lyapunov and Sylvester Solver

debug = false;

A = varargin{1};
B = varargin{2};
C = varargin{3};

if nargin > 3
    poles = varargin{4};
end

if ~isinf(max(size(A)))
    % error('cqtlyap is only supported for infinite matrices');
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
    [xm, xp] = evinterp(@(a,b,c) c ./ (a + b), am, ap, bm, bp, -cm, -cp);
end

% Solve the correction equation
X = cqt(xm, xp, [], [], size(A, 2), size(B, 1));
R = A*X + X*B + C;

[RU, RV] = correction(R);

% Solve the equation for the correction
if ~exist('poles', 'var')
    [Xu, Xv] = ek_sylv(A, B, RU, RV, inf, ...
        cqtoption('threshold'), debug, inf);
else
    [Xu, Xv] = rk_sylv(poles, A, B, RU, RV, inf, ...
        cqtoption('threshold'), debug, inf);
end

X = X + cqt([], [], Xu, Xv);

end

