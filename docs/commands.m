%% Commands available in CQT
%
% This page describes all commands that are available in the CQT Toolbox;
% documentation for each of these commands is available by clicking on
% their name. 
%

%% Arithmetic operations
%
% * |C = A + B| and |C = A - B| compute the sum and differences of QT
%   matrices.
% * |C = A*B| computes the matrix product between |A| and |B|. 
% * |C = A \ B| computes the action of the inverse of |A| on |B|. Note that
% |B| may be a vector obtained with |B = cqt(v)| (represented as a QT
% matrix with zero symbol), which allows to solve linear systems with QT
% matrices. See <LinearSystems.html LinearSystems> for further details. 

%% General commands
%
% * <cmd_cond.html cond>: The |cond| command computes (or estimates) the condition number 
%    of a QT matrix.
% * <cmd_conj.html conj>: Returns the entry-wise complex conjugate of a QT
%    matrix.
% * <cmd_correction.html correction>: Extracts the correction part of a QT
%    matrix. 
% * <cmd_cqtinfo.html cqtinfo>: Prints brief information on a QT matrix. 
% * <cmd_cqt.html cqt>: Constructs a new QT matrix. 
% * <cmd_cqtoption.html cqtoption> sets global options for the QT toolbox. 
% * <cmd_cqtrank.html cqtrank>: Returns the rank of the compact correction.
% * <cmd_ctranspose.html ctranspose>: Returns the
%   complex-conjugate-transpose of the QT matrix $A$. 
% * <cmd_diag.html diag> returns the a QT matrix that only keeps diagonal
%   entries from the original one. 
% * <cmd_eye.html eye> constructs an identity matrix in QT format. 
% * <cmd_full.html full> returns the full dense version of a QT matrix of
%   finite size. 
% * <cmd_inv.html inv> returns the inverse of a QT matrix, again in the QT
%   format. 
% * <cmd_limit.html limit> computes a vector with the limits of every
%   column for large row indices. 
% * <cmd_normest.html normest> Norm estimator, based on a power method. 
% * <cmd_norm.html norm> Computation of the norm of QT matrices. 
% * <cmd_size.html size> returns the size of a QT matrix. 
% * <cmd_symbol.html symbol> accesses the symbol of the QT matrix under
%   consideration. 
% * <cmd_transpose.html transpose> computes the transpose of a QT matrix
% * <cmd_ul.html ul> Computes the UL factorization, that is used (among 
%   other things) to solve linear systems and compute the inverse. 

%% Eigenvalue computations
%
% * |basins|: Compute the basin of attraction for the fixed point iteration
%   used when computing eigenvalues of QT matrices. 
% * |distances|:
% * |eig_all|: 
% * |eig_single|
% * |range|

%% Matrix equations and matrix functions
%
% * |cqtlyap|
% * |cqtstein|
% * |expm|
% * |funm|
% * |logm|
% * |polyvalm|
% * |sqrtm|