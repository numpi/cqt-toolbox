function c = cond(A, norm_type)
%COND Compute (or estimate) the condition number of A
%
%     C = COND(A) computes the condition number of the CQT
%     matrix A. The norm used is the 2-norm for finite matrices,
%     and the CQT-norm for semiinfinite ones.
%
%    C = COND(A, NORM_TYPE) computes the condition number for the
%    given norm. Valid choices for finite matrices are NORM_TYPE
%    set to 1, inf, 2; for semiinfinite matrices 'CQT' is the only
%    valid choice.

if ~exist('norm_type', 'var')
	c = norm(A) * norm(inv(A));
else
	c = norm(A, norm_type) * norm(inv(A), norm_type);
end
