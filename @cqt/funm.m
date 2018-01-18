function F = funm(A, fun, varargin)
%FUNM computes the matrix function FUN on the argument M
%
% F = FUNM(A, FUN) computes the function FUN on the input matrix argument A
%
% F = FUNM(A, FUN, 'opt1', value1, 'opt2', value2, ...) performs the 
% computation with the optional parameters opt1, opts2, ... with the given
% values. Possible choices for the parameters are: 
%
% - 'max_it': Maximum number of points for the contour integral formula.
%             The default value is 64. 
% - 'poles':  Poles of the functions, if any. Should be specified as a
%             vector, and needs to be given for the formula to work if any
%             of them is smaller than the norm of A. 

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
if max(A.sz) == inf
        r = norm(A) + 1;
else
        r = norm(A, inf) + 1;
end

while sum(abs(poles) == r) > 0
	r = r + 1;
end

F = contour_integral(res, 0, norm(A)+1, max_it);

for i = 1 : length(poles)
	% If the pole is inside the radius, remove the integral around it
	if abs(poles(i)) < r
		corr = inv(cqt(poles(i), poles(i), [], [], m, n) - A);
		rr = .5 / norm(corr, inf);
		
		corr = contour_integral(res, poles(i), rr, max_it);
		F = F - corr;
	end
end

