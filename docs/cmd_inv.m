%% Inverse of a QT matrix
%
% The command |inv| returns the inverse of a QT matrix, represented in the
% QT format. To compute the inverse, it is necessary that the symbol $a(z)$
% defining the Toeplitz part needs to be invertible, i.e., it needs to
% never vanish over the unit circle. 

%% Syntax
% * |iA = inv(A)|

%% Example
%
% We choose a well-conditioned Toeplitz matrix, that can be inverted
% easily. 

A = cqt([ 3 1 ], [ 3  1 ]);
iA = inv(A);
cqtinfo(iA)

%%
% We can then test the accuracy of the computed inverse by testing that the
% inverse times the original matrix is close to the identity; we expect
% this number to be close to the truncatio threshold, that can be set with
% the <cmd_cqtoption.html cqtoption> command. 

res = norm( A * iA - eye(size(A), 'like', cqt), 'cqt')