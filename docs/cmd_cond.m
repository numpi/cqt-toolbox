%% Condition number
%
% The |cond| command computes the condition number
% $\| A \| \cdot \| A^{-1} \|$ of a QT matrix.
%
% The default is to compute 
% the condition number with respect to the spectral norm for finite 
% matrices, and with respect to the QT norm for infinite ones. 
%
% Another choice of norm can be given using the syntax |cond(A, nrm_type)|,
% where |nrm_type| is some norm supported by the <norm.html norm> command. 

%% Syntax
% * |c = cond(A)|
% * |c = cond(A, nrm_type)|

%% Example
%
% The Toeplitz matrix with symbol $a(z) = z^{-1} + 3 + z$ is always well
% conditioned, since $a(z)$ is far from zero on the unit circle. 

A = cqt([3 1], [3 1]);
cond(A)