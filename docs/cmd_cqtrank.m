%% Rank of the correction
%
% The command |cqtrank| returns the rank of the top-left correction in the
% QT matrix; for matrices of finite size, the rank of the bottom right
% correction is returned as a second output. 

%% Syntax
% * |K = cqtrank(C)|
% * |[KT, KB] = cqtrank(C)|

%% Example
%
% The rank of the correction in the inverse of a tridiagonal Toeplitz
% matrices is at most one. 

A = inv(cqt([3 1], [3 1]));
rk = cqtrank(A)