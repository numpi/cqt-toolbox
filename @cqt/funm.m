function F = funm(A, fun, varargin)
% FUNM computes the matrix function FUN on the argument M
%
% F = FUNM(A, FUN) computes the function FUN on the input matrix argument M
%

p = inputParser;

addParamValue(p, 'max_it', 64, @isscalar);
addParamValue(p, 'poles', [], @isvector);

parse(p, varargin{:});

max_it = p.Results.max_it;
poles = p.Results.poles;

[m, n] = size(A);

if m ~= n
	error('The matrix A needs to be square');
end

res = @(z) inv(cqt(z, z, [], [], m, n) - A) * fun(z);

% Make sure that we do not choose a circle which contains poles for the
% contour integral.
r = norm(A, inf) + 1;

while sum(abs(poles) == r) > 0
	r = r + 1;
end

F = contour_integral(res, 0, norm(A)+1, max_it);

for i = 1 : length(poles)
	% If the pole is inside the radius, remove the integral around it
	if abs(poles(i)) < r
		corr = inv(cqt(poles(i), poles(i), [], [], m, n) - A);
		rr = .5 / norm(corr, inf)
		
		corr = contour_integral(res, poles(i), rr, max_it);
		F = F - corr;
	end
end

