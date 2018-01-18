function Y = polyvalm(p,X)
%POLYVALM Evaluate polynomial with matrix argument.
%
%   Y = POLYVALM(P,X), when P is a vector of length N+1 whose elements
%   are the coefficients of a polynomial, is the value of the
%   polynomial evaluated with matrix argument X.  X must be a
%   square matrix.
%
%       Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)*I
%
% Polynomial evaluation p(x) using Horner's method.

% Check input is a vector
if ~(isvector(p) || isempty(p))
	error('POLYVALM: first argument must be a non empty vector');
end

np = length(p);
[m, n] = size(X);

if m ~= n
	error(message('Polyvalm: NonSquareMatrix'))
end

if np == 1    %Quick return if possible.
	Y = cqt(p(1), p(1), [], [], m, n);
	return
end

Y = cqt([],[],[],[],m,n);
for i = 1:np
	Y = X * Y + cqt(p(i),p(i), [], [], m, n);
end

