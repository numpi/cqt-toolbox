%% |cqtinfo|
%
% The |cqtinfo| commands prints brief information on a QT object, such as
% the length of the symbol, and the support and rank of the correction. 

%% Syntax
% * |cqtinfo(A)|

%% Example
%
% We test the command on a infinite Toeplitz matrix with a random 
% correction of rank $2$. 

A = cqt([1 2], 1, rand(12, 2), rand(8, 2));
cqtinfo(A)
