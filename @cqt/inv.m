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

        % The above correction may given an output with two or more rows
        % due to rounding effects; since the above should be a row vector,
        % we ignore any other data
        if ~isempty(ciT)
            ciT = ciT(1, :);
        end
        
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
