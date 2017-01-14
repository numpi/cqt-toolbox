function [l, u] = reciprocal(vm, vp, nm, np)
%RECIPROCAL Compute the inverse of a Laurent polynomial. 

if ~exist('nm', 'var')
	nm = inf;
end

if ~exist('np', 'var')
	np = inf;
end

if length(vm) == 1
	u = reciprocal_poly(vp, np);
	l = u(1);
elseif length(vp) == 1
	l = reciprocal_poly(vm, nm);
	u = l(1);
else
	switch cqtoption('inversion')
		case 'cr'
			[l, u] = reciprocal_cr(vm, vp, nm, np);
		case 'fft'
			[l, u] = reciprocal_fft(vm, vp, nm, vp);
	end
end



end

