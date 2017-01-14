function q = reciprocal_poly(p, d)
%RECIPROCAL Compute the Taylor expansion of the inverse of a polynomial.
%
% Q = RECIPROCAL_POLY(P, D) computes at most N entries of the Taylor series
% of the inverse of the polynomial 
%
%   P(1) + P(2) * x + ... + P(l+1) x^l
%

if isempty(p) || p(1) == 0
	error('The given polynomial is not invertible as series');
end

if length(p) == 1
	q = 1 / p;
	return;
end

n = length(p) - 1;

T = toep(p(1), p(1:end-1), n, n);
A = toep(p(end:-1:2), p(end), n, n);

q = [ 1, zeros(1, n-1) ] / T;
accurate = false;

threshold = norm(q, inf) * eps;

while length(q) < d && ~accurate
	nq = -q(end-n+1:end) * A / T;
	
	if norm(nq, inf) < threshold
		accurate = true;
	end
	
	q = [ q , nq ];
end

q = truncate(cln(q), d);

end

