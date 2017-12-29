function [G, R, B0] = cr(Am1, A0, A1, max_it)
%CR	Executes the cyclic reduction iteration for solving Am1 + A0 X + A1 X^2 = 0
%
%	[G, R] = CR(Am1, A0, A1) computes the matrices G and R solutions of
%	Am1 + A0 X +A1 X^2 = 0 and A1 + A0 X +Am1 X^2 = 0, respectively
%
%	[G, R] = CR(Am1, A0, A1, MAX_IT) limits the number of iterations to MAX_IT,
%	the default value is MAX_IT=20
%
%	[G, R, B0] = CR(Am1, A0, A1) returns G, R and the limit point of the sequence A0^(h)
%
%	[G, R, B0] = CR(Am1, A0, A1, MAX_IT) limits the number of iterations to MAX_IT,
%	the default value is MAX_IT=20

if ~exist('max_it','var')
	max_it = 20;
end

Bm1 = Am1;
B0 = A0;
B1 = A1;

hB0 = B0;

if size(A0,1) == inf
	norm_type ='cqt';
else
	norm_type = inf;
end

tol = 1e-8;

for i=1:max_it
	BB0 = inv(B0);
	temp1 = BB0 * B1;	
	temp2 = BB0 * Bm1;
	temp3 = Bm1 * temp1;
	temp4 = B1 * temp2;
	hB0 = hB0 - temp4;
	B0 = B0 - temp4 - temp3;
	Bm1 = -Bm1 * temp2;
	B1 = -B1 * temp1;

	% norm(Bm1, norm_type)
	% norm(B1, norm_type)
	% norm(symbol(Bm1))
	% s = svd(correction(Bm1)); 
	% s(1:min(3,length(s)))
    % min(norm(Bm1, norm_type), norm(B1, norm_type))
	% G = hB0 \ Am1;
	% S = Am1  + A0 * G + A1
	% max(max(abs(correction(B1))))
    
%    [BU,BV] = correction(B0);
%    if size(BU,1) > 20 * size(BV,1)
%        B0 = extrapolate_limit(B0, tol);
%        B1 = extrapolate_limit(B1, tol);
%        Bm1 = extrapolate_limit(Bm1, tol);
%        hB0 = extrapolate_limit(hB0, tol);
%    end
    
%    fprintf('Iteration %d, res = %e', i, min(norm(Bm1, norm_type), norm(B1, norm_type)));
    
	if min(norm(Bm1, norm_type), norm(B1, norm_type)) < eps
		break
	end
end
if i==max_it
	warning('maximum number of iterations reached!')
end

G = - hB0 \ Am1;
R = - A1 / hB0;
