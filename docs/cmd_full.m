%% Converting finite QT matrices to dense matrices
%
% The command |full| returns a dense representation of a QT matrix,
% similarly to what the MATLAB command does for sparse matrices. Clearly,
% this only makes sense for matrices of finite size. 

%% Syntax
% * |F = full(A)|

%% Example
%

A = cqt([1 2 3], [1 4], rand(2), rand(3), 7, 7);
full(A)