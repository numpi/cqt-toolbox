function [ iT ] = inv(T)
%INV Inverse of a CQT-matrix T.
%
% iT = INV(T) computes the inverse of a CQT matrix T and produces a new CQT
% matrix iT.
%

if T.sz(2) == inf && T.sz(1) == T.sz(2)
    [n, p, U, V] = qt_inv(T.n, T.p, T.U, T.V);
    iT = cqt(n, p, U, V);
    
    if ~isempty(T.c)
        ciT = correction(cqt(T.c) * iT);
        
        S = extend(cqt([]), ciT);
        S = iT * S / (1 + sum(ciT));
        iT = iT - S;
    end
elseif T.sz(1) == T.sz(2)
    [n, p, U, V, W, Z] = fqt_inv2(T.n, T.p, T.U, T.V, T.W, T.Z, T.sz(1));
    iT = cqt(n, p, U, V, W(end:-1:1,end:-1:1), ...
        Z(end:-1:1,end:-1:1), T.sz(1), T.sz(2));
else
    error('Can not invert a non square matrix');
end
end
