function [ T ] = mtimes(T1, T2)
%MTIMES Multiply two CQT matrices T1 and T2. 
%
% T = MTIMES(T1, T2) computes the CQT matrix T obtained multiplying two CQT
%     matrices T1 and T2. If T2 is a dense finite matrix it is considered
%     as a semiinfinite matrix with only the leading top-left corner
%     different from zero. 

if isa(T2, 'cqt')
    [cm,cp,cU,cV]=qt_mult(T1.n, T1.p, T1.U, T1.V, T2.n, T2.p, T2.U, T2.V);
    T = cqt(cm, cp, cU, cV);
else
    [cm,cp,cU,cV]=qt_mult(T1.n, T1.p, T1.U, T1.V, ...
        0, 0, T2, eye(size(T2,2)));
    T = cqt(cm, cp, cU, cV);
end


end

