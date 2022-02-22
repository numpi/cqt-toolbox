%% Solving linear systems with QT matrices
%
% Given a QT matrix $A$, we may want to solve a linear system $Ax = b$.
% Since $A$ is infinite-dimensional, the same holds for $x$ and $b$, and we
% need therefore to assume that they can be truncated to a finite length. 
%
% If $b$ has finite support, we can represent it as a vector of infinite
% length by slightly "abusing" the QT representation, and interpret it as a
% QT matrix with zero symbol, and a correction in the first column only. 

b = cqt([ 1 ; 2 ; 3 ; 4 ]);

%%
% We now define a matrix $A$. If the symbol is invertible, we can evaluate
% $x = A^{-1}b$, and this is overloaded in MATLAB by means of the backslash
% operator. 

A = cqt([ 4 1 ], [ 4 -1], rand(4)); % QT matrix with random correction
x = A \ b;

%%
% Right now |x| is a QT matrix as well, again with zero-symbol and a
% correction in the first column. To extract the solution of the linear
% system truncated to a finite length vector, we use the |correction|
% command:

fx = correction(x); length(fx)

%% 
% We can check the norm of the residual in the 2-norm by using the norm 
% command; this should be comparable to the current accuracy set with
% |cqtoption('threshold')|. Note that all these object are of infinite
% size, and truncation is performed automatically. 

norm(A * x - b, 2) / norm(b, 2)

%% 
% In practice, one can also check that the residual is small on a "large"
% enough finite section of $A$. If $n$ is not chosen large enough, then the
% result will be inaccurate. 

n = 10;
norm(A(1:n, 1:n) * x(1:n, 1) - b(1:n, 1)) / norm(b(1:n, 1))

n = 1e3;
norm(A(1:n, 1:n) * x(1:n, 1) - b(1:n, 1)) / norm(b(1:n, 1))

%% UL factorizations
%
% In practice, linear systems are not solved by explicitly computing the
% inverse. Rather, we factorize $A$ into an UL factorization plus a
% low-rank correction, by computing matrices U, L, E such that 
%
% $$ A = UL + E $$
%
% where $U$ is upper triangular and Toeplitz, $L$ is lower triangular and
% Toeplitz, and $E = WZ^*$ is a low-rank correction. Without the
% correction, this allows to solve linear systems by a (properly truncated)
% backsubstitution scheme. When $E \neq 0$, we use a Sherman-Morrison-Woodbury 
% formulate to evaluate $x = A^{-1}b$. 
%
% The factorization can be computed in MATLAB by using the |ul| command:

[U, L, E] = ul(A);
norm(U*L + E - A)