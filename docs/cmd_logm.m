%% Matrix logarithm of QT matrices
%
% The function |logm| evaluates the matrix logarithm using inverse scaling
% and squaring, and a Padè approximant. 

%% Syntax
% |X = logm(A)|

%% Example
% 
% For the logarithm to be well-defined, we need to avoid negative real
% eigenvalues; here, we choose a positive definite matrix, and check that
% the matrix exponential of the outcome of |logm| is again the original
% matrix. 

A = cqt([ 3 1 ], [ 3 1 ]);
L = logm(A);
res = norm(expm(L) - A) / norm(A)