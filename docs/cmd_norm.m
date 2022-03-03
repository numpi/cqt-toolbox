%% Computation of norms of QT matrices
%
% The |norm| command can be used to evaluate the norm of a QT matrix.
% Several norms are available, and can be selected with the second optional
% argument to the command. The default is to use the QT norm for infinite
% matrices, and the spectral norm for finite ones. 
%
% The implemented norms, selected by the |nrm_type| parameter, should be
% one of |'cqt'|, |'eqt'| (for extended matrices) or |p|, where |p = 1, 2, inf|. 

%% Syntax
% * |nrm = norm(A)|
% * |nrm = norm(A, nrm_type)|

%% Example
%
A = cqt([ 3 1 ], [ 3 2 ])

nrmCQT = norm(A, 'cqt')
nrmInf = norm(A, inf)
nrm1 = norm(A, 1)
