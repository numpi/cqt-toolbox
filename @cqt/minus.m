function T = minus(T1, T2)
%MINUS subtract two CQT matrices.
%
%     T = MINUS(T1, T2) subtracts the CQT-matrix T2 from  the CQT-matrix T1
%     and produces a new CQT matrix T.

if isa(T1, 'cqt') && isa(T2, 'cqt')
    T = T1 + (-T2);
elseif ~isa(T1, 'cqt')
    error('Incompatible types addition. \nIf you wan to subtract ' + ... 
        'a matrix of %s with a cqt matrix T you can' + ...
        'use cqt(A) - T', class(T1));
elseif ~isa(T2, 'cqt')
    error('Incompatible types addition. \nIf you wan to subtract ' + ...
        'a cqt matrix T  with a matrix of %s you can' + ...
        'use T - cqt(A)', class(T2));
end
end
