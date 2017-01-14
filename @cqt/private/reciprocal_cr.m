function [l, u] = reciprocal_cr(vm, vp, nm, np)
%RECIPROCAL_CR Compute the inverse power series of VM and VP. if NM and NP
%are specified, at most NM and NP coefficients are computed.

if ~exist('nm', 'var')
	nm = inf;
end

if ~exist('np', 'var')
	np = inf;
end

accurate = false;

kleft = 1;
kright = 1;

vm = reshape(vm, 1, length(vm));
vp = reshape(vp, 1, length(vp));

threshold = eps * ( norm(vm) + norm(vp) );

while ~accurate && ( kleft <= nm || kright <= np )
	% Increase the padding to gain more accuracy
	kleft  = 2 * kleft;
	kright = 2 * kright;
	
	[~,~,l,u] = spectral_cr([ vm , zeros(1, kleft) ], ...
		[ vp, zeros(1, kright) ]);
	
	if abs(l(end)) < threshold && abs(u(end)) < threshold
		accurate = true;
	end
end

l = truncate(cln(l), nm).';
u = truncate(cln(u), np).';

end

