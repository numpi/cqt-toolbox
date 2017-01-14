function [n, p] = symbol(T)
%SYMBOL Returns the symbol of the Toeplitz part in the CQT matrix T.
%
%   [N, P] = SYMBOL(T) returns two vectors containing the negative and
%       positive coefficients of the symbol of the Toeplitz part of T,
%       respectively.

n = T.n;
p = T.p;

end

