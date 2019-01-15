function X = cqtstein(A, B, C, poles)
%CQTSTEIN Solves the equation A*X*B + X = C
%
% X = CQTSTEIN(A, B, C) computes the CQT solution to the matrix equation 
%
%    (1)   A*X*B + X = C 
%
%     under the assumption that both A and B have spectral radius strictly
%     smaller than 1. 
%
% X = CQTSTEIN(A, B, C, POLES) specifies the poles to use in the Sylvester
% equation constructed by left and right multiplying by shifted inverses of
% A and B. 

M1 = inv(A + cqt(1,1));
M2 = inv(cqt(1,1) - B);

if ~exist('poles', 'var')
    X = cqtlyap(M1, B * M2, - M1 * C * M2);
else
    X = cqtlyap(M1, B * M2, - M1 * C * M2, poles);
end

end