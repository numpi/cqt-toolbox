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
% * <cmd_ul.html ul> computes the UL factorization, that is used (among 
%   other things) to solve linear systems and compute the inverse. 

%% Eigenvalue computations
%
% * <cmd_basins.html basins> computes the basin of attraction for the
%   matrix iterations used when computing eigenvalues of QT matrices. 
% * <cmd_distances.html distances> computes the distance between the
%   eigenvalues of finite truncations of the infinite-dimensional operator,
%   and its true eigenvalues. 
% * <cmd_eig_all.html eig_all> computes the spectrum of a QT matrix
%   (decomposed as the continuous part and the isolated eigenvalues). 
% * <cmd_eig_single.html eig_single> computes a single isolated eigenvalue
%   of a QT matrix starting from an initial approximation. 
% * <cmd_range.html range> plots the boundary of the continuous spectrum of
%   $A$, which coincides with the symbol evaluated at all the points in the
%   unit circle. 

%% Matrix equations and matrix functions
%
% * <cmd_cqtlyap.html cqtlyap> solves Lyapunov equations $AX + XA + C = 0$
%   and Sylvester equations $AX + XB + C = 0$. 
% * <cmd_cqtstein.html cqtstein> solves the Stein equation $AXA + X = C$. 
% * <cmd_expm.html expm> computes the matrix exponential of $A$. 
% * <cmd_funm.html funm> evaluates a generic matrix function by means of a
%   contour integral. 
% * <cmd_logm.html logm> evaluates the matrix logarithm of $A$. 
% * <cmd_polyvalm.m polyvalm> evaluates a polynomial at a matrix $A$. 
% * <cmd_sqrtm.m sqrtm> evaluate the matrix square root of $A$. 