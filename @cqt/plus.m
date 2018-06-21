function T = plus(T1, T2)
%PLUS Add two CQT matrices.
%
% T = PLUS(T1, T2) adds two CQT matrices and produces a new CQT matrix T.

if isa(T1, 'cqt') && isa(T2, 'cqt')
    if T1.sz(1) ~= T2.sz(1) || T1.sz(2) ~= T2.sz(2) 
        error('Incompatible dimensions in the addition')
    end
    if max(T1.sz) == inf
        [ cm, cp, cu, cv ] = qt_add(T1.n, T1.p, T1.U, T1.V, ...
            T2.n, T2.p, T2.U, T2.V);
        T = cqt(cm, cp, cu, cv, [], [], T1.sz(1), T1.sz(2));
    else
        [ cm, cp, cu, cv, cw, cz ] = fqt_add2(T1.n, T1.p,  ...
            T1.U, T1.V, T1.W, T1.Z, T2.n, T2.p, T2.U, ...
            T2.V, T2.W, T2.Z);
        T = cqt(cm, cp, cu, cv);
        T.W = cw; T.Z = cz; T.sz = T1.sz;
        T = merge_corrections(T);
    end
    
    T.c = formatted_sum(T1.c, T2.c);
elseif ~isa(T1, 'cqt')
    error('Incompatible types addition. \nIf you want to add a' + ...
        'matrix of %s with a cqt matrix T you can use cqt(A) + T', ...
        class(T1));
elseif ~isa(T2, 'cqt')
    error('Incompatible types addition. \nIf you want to add a' + ...
        'cqt matrix T  with a matrix of %s you can use T + cqt(A)', ...
        class(T2));
end
end

