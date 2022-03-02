%% Complex conjugate
%
% The |conj| command returns the complex conjugate of a QT matrix $A$; no
% transposition is performed, and the operation is carried out entry-wise.
% For the complex conjugate transpose, use |A'| instead. 

%% Syntax
% * |Ac = conj(A)|

%% Example
%
A = cqt([ 1 + 1i, 2 + 1i ], 1+ 1i);
B = conj(A); 
B(1:4, 1:4)