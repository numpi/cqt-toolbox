function [ iT ] = inv(T)
%INV Inverse of a CQT-matrix T. 
%
%     iT = INV(T) computes the inverse of a CQT matrix T and produces a new CQT matrix iT.
%
[n, p, U, V] = qt_inv(T.n, T.p, T.U, T.V);
iT = cqt(n, p, U, V);
