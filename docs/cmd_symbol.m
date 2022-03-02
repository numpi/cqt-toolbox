%% Symbol of a QT matrix
%
% The command |symbol| returns the coefficients of the Laurent polynomial
% that approximates the symbol of a QT matrix. 
%
% The correction part can be obtained with the <cmd_correction.html
% correction> command. 

%% Syntax
% * |[am, ap] = symbol(A)|

%% Example
%
% The two vectors returned by the |symbol| command contain the negative and
% positive coefficients of the Laurent polynomial for the Toeplitz part of
% $A$. 

A = cqt([ 1 2 3 ], [ 1 4 5 ]);
[am, ap] = symbol(A)