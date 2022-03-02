%% Norm estimation
%
% The matrix-vector product by a QT matrix $A$ can be computed cheaply, and
% therefore the norm of $QT$ matrices is amenable to be computed by a power
% iteration. 
%
% The estimate norm is the spectral norm, and should be used to get rough
% estimates. Note that in general <cmd_norm.html norm> is relatively cheap
% for QT matrices, and can be a much more accurate alternative. 

%% Syntax
% * |nrm = normest(A)|

%% Example
%
% We use normest, and compare the result with the "exact" norm compute by
% the <cmd_norm.html norm> command. 

A = cqt([3 2 1], [3 1]);
nrmA = norm(A, 2)
nrmEst = normest(A)