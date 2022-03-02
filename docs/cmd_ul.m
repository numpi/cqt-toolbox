%% UL Factorization
%
% QT matrices can be factorized in the form $A = UL + E$, where $U$ is
% Toeplitz and upper triangular, $L$ is Toeplitz and lower triangular, and
% $E$ is a compact correction. The function |ul| returns this
% factorization. 

%% Syntax
% * |[U,L,E] = ul(A)|

%% Example
%

A = cqt([ 3 1 ], [ 3 2 ], rand(4));
[U, L, E] = ul(A);
res = norm(U*L + E - A) / norm(A)