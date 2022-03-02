%% Size of a QT matrix
%
% The command |size| returns the number of rows and columns of a QT matrix.
% Note that these can be infinity. 

%% Syntax
% * |sz = size(A)|
% * |sz = size(A, dim)|

%% Example
%

A = cqt(1, 1); 
szA = size(A)