function T = plus(T1, T2)
%PLUS Add two CQT matrices. 
%
% T = PLUS(T1, T2) adds two CQT matrices and produces a new CQT matrix T. 

if isa(T2, 'cqt')
    [ cm, cp, cu, cv ] = qt_add(T1.n, T1.p, T1.U, T1.V, ...
        T2.n, T2.p, T2.U, T2.V);
    T = cqt(cm, cp, cu, cv);
else
    [U,S,V] = svd(T2);
    [ cm, cp, cu, cv ] = qt_add(T1.n, T1.p, T1.U, T1.V, ...
        0, 0, U * S, conj(V));
    T = cqt(cm, cp, cu, cv);
end


end

