%% Tranpose of a QT matrix
%
% The |transpose| command, which can be accessed using the syntax |A.'|,
% constructs a representation of the transpose of the matrix |A| (without
% any complex conjugation). 

%% Syntax
% * |B = A.'|
% * |B = transpose(A)|

%% Example

A = cqt([1, 2+1i], 1);
A.'