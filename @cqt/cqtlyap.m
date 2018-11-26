function X = cqtlyap(A, B, C)
%CQTLYAP Lyapunov and Sylvester Solver

debug = false;

if ~isinf(max(size(A)))
    error('cqtlyap is only supported for infinite matrices');
end

% Solve the equation on the Toeplitz partby evaluation / interpolation
[am, ap] = symbol(A);
[bm, bp] = symbol(B);
[cm, cp] = symbol(C);

% Construct the symbol of X
[xm, xp] = evinterp(@(a,b,c) c ./ (a + b), am, ap, bm, bp, -cm, -cp);

% Solve the correction equation
X = cqt(xm, xp);
R = A*X + X*B + C;

[RU, RV] = correction(R);

% Solve the equation for the correction
[Xu, Xv] = ek_sylv(A, B, RU, RV, inf, ...
    cqtoption('threshold'), debug, 'cqt');

X = X + cqt([], [], Xu, Xv);

end

