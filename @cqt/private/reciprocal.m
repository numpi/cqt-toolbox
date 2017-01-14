function [l, u] = reciprocal(vm, vp, nm, np)
%RECIPROCAL Compute the inverse of a Laurent polynomial. 
%
% [L, U] = RECIPROCAL(VM, VP) computes the coefficients of the Laurent
% polynomial with positive coefficients VP and negative coefficients VM.
% The coefficients are order starting from the constant to the one of
% higher (negative or positive) degree. This function computes all the
% non-negligible coefficients with respect to machine precision. 
%
% [L, U] = RECIPROCAL(VM, VP, NM, NP) does the same computation of
% RECIPROCAL(VM, VP) but only returns at at most NM negative coefficients
% and NP positive ones. 
%
% Authors: Stefano Massei <stefano.massei@sns.it>
%          Leonardo Robol <leonardo.robol@cs.kuleuven.be>

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

