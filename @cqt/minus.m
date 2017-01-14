function T = minus(T1, T2)
%MINUS subtract two CQT matrices. 
%
%     T = MINUS(T1, T2) subtracts the CQT-matrix T2 from  the CQT-matrix T1 
%     and produces a new CQT matrix T. 

if isa(T1, 'cqt') && isa(T2, 'cqt')
    if max(T1.sz) == inf
        [ cm, cp, cu, cv ] = qt_add(T1.n, T1.p, T1.U, T1.V, ...
            -T2.n, -T2.p, -T2.U, T2.V);
        T = cqt(cm, cp, cu, cv, [], [], T1.sz(1), T1.sz(2));
    else
        [ cm, cp, cu, cv, cw, cz ] = fqt_add2(T1.n, T1.p,  ...
        		T1.U, T1.V, T1.W, T1.Z, -T2.n, -T2.p, -T2.U, ...
			T2.V, -T2.W, T2.Z);
    		T = cqt(cm, cp, cu, cv);
 		T.W = cw; T.Z = cz; T.sz = T1.sz;
    end
elseif ~isa(T1, 'cqt')
    error('Incompatible types addition. \nIf you wan to subtract a matrix of %s with a cqt matrix T you can use cqt(A) - T',class(T1));
elseif ~isa(T2, 'cqt')
    error('Incompatible types addition. \nIf you wan to subtract a cqt matrix T  with a matrix of %s you can use T - cqt(A)',class(T2));
end
end
