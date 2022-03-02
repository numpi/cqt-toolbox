%% Correction extraction
%
% QT matrices are defined as $A = T(a(z)) + E_a$, where $E_a$ is a compact
% correction approximated by a low-rank matrix with finite support. 
%
% This correction can be extracted, either as a full matrix or in factored
% form (the latter is usually much better from the efficiency viewpoint). 
%
% The symbol can be accessed using the <cmd_symbol.html symbol> command. 

%% Syntax
% * |E = correction(A)| 
% * |[U, V] = correction(A)|

%% Example
%
% A simple way to generate matrices with a correction is taking the inverse
% of a Toeplitz matrix which is not triangular. 

A = cqt([3 1], [3 1]);
iA = inv(A);

%%
% The correction can be extracted in full form

EA = correction(iA); 
support = size(EA)
EA(1:4, 1:4)

%%
% ... or in factored form, which requires much less storage (notice how in
% this case the correction has rank equal to $1$). 
[UA, VA] = correction(iA);
supportU = size(UA)
supportV = size(VA)