function X = cqtstein(A, B, C, poles)
%CQTSTEIN Solves the equation A*X*B + X = C
%
% X = CQTSTEIN(A, B, C) computes the CQT solution to the matrix equation 
%
%    (1)   A*X*B + X = C 
%
%     under the assumption that both A and B have spectral radius strictly
%     smaller than 1. 

I = cqt(1,1);

M1 = inv(A - I);
M2 = inv(B + I);

X = cqtlyap(M1 * (A + I), M2 * (B - I), 2 * M1 * C * M2);

end