function r = norm(T)
%NORM CQT-norm of a CQT-matrix
%
%     r = NORM(T) computes the CQT-norm of a CQT-matrix  

r = qt_norm(T.n, T.p, T.U, T.V);

