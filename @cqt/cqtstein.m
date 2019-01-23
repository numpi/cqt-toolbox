function X = cqtstein(A, B, C, varargin)
%CQTSTEIN Solves the equation A*X*B + X = C
%
% X = CQTSTEIN(A, B, C) computes the CQT solution to the matrix equation 
%
%    (1)   A*X*B + X + C = 0
%
%     under the assumption that both A and B have spectral radius strictly
%     smaller than 1. 

p = inputParser;

addParameter(p, 'tol', cqtoption('threshold'));

parse(p, varargin{:});

tol = p.Results.tol;

IA = eye(size(A), 'like', A);
IB = eye(size(B), 'like', B);

M1 = inv(A - IA);
M2 = inv(B + IB);

X = cqtlyap(M1 * (A + IA), M2 * (B - IB), 2 * M1 * C * M2, 'tol', tol);

end