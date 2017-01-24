function X = fqt_solvetriang(A, B)
%FQT_BACKSUB Compute A \ B, where A is triangular and B a dense matrix.

n = size(B, 1);

if min(length(A.n), length(A.p)) ~= 1
	error('The matrix A is not triangular');
end

if length(A.n) == 1
	sA = spdiags(ones(n, 1) * A.p, 0:(length(A.p)-1), n, n);
else
	sA = spdiags(ones(n, 1) * A.n, 0:-1:-(length(A.n)-1), n, n);
end

X = sA \ B;

end

