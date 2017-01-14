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
	[l, u] = spectral(vm, vp, nm, np);
	
	[linv, ~] = reciprocal_poly(l, nm);
	[~, uinv] = reciprocal_poly(u, np);
	
	r = conv_fft(linv, [ zeros(1, length(linv) - 1), uinv ]);
	l = r(1:length(linv));
	u = r(length(linv):end);
end



end

