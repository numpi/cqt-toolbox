%% Arithmetic with infinite matrices
%
%% The first steps with QT matrices
% The CQT toolbox allows to represent infinite Quasi-Toeplitz matrices
% directly in MATLAB. Assume that we have a semi-infinite Toeplitz matrix
% $T$ with symbol $s(z)$ defined as follows:
%
% $$ s(z) = 2 z^{-2} + z^{-1} - 3 - z + 4z^2 $$
%
% We can represent this matrix using the CQT toolbox as

T = cqt([ -3 1 2 ], [ -3 -1 4 ])

%% Adding a finite perturbation
% We can add a perturbation in the top-left corner, by specifying it as a
% dense $m \times n$ matrix E.

E = randn(3, 4);
T = cqt([ -3 1 2 ], [ -3 -1 4 ], E)

%% Accessing the elements of the matrix
%
% The matrix we have defined is of infinite size. We can however access
% finite sections of it with the familiar notation used in MATLAB to slice
% arrays and matrices:

T_top = T(1:6, 1:6)
