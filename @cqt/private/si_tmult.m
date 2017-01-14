function [cm,cp,cU,cV]=si_tmult(am, ap, bm, bp)
%SI_TMULT Compute the product of two semi-infinite Toeplitz matrices.
%
% Compute the product C=AB of two semi-infinite Toeplitz matrices
% A and B defined by their first column and first row am,ap
% and bm,bp, respectively
% the matrix product C is written as toep(cm,cp) + cU*cV

nam = length(am); nap = length(ap);
nbm = length(bm); nbp = length(bp);

if nam == 0 || nbm == 0
	cm = [];
	cp = [];
	cU = [];
	cV = [];
	return;
end

% resize the input vectors according to their numerical length
[nam_eps,am] = epslength(am);
[nap_eps,ap] = epslength(ap);
[nbm_eps,bm] = epslength(bm);
[nbp_eps,bp] = epslength(bp);

% Convolution c(z) of a(z) e b(z), with detection of the center
a = [am(end:-1:1),ap(2:end)]; centera = length(am);
b = [bm(end:-1:1),bp(2:end)]; centerb = length(bm);
c = conv(a,b); m = centera + centerb-1;
cm = c(m:-1:1);cp = c(m:end);
% clean cm and cp
cm = cln(cm);
cp = cln(cp);
% Compute R
na = length(am); nb = length(bp);
if na <=1 || nb <= 1
	cU = []; cV = [];
else
	if na>nb
		am = am(2:end);   bm = [bp(2:end),zeros(1,na-nb)];  bp = zeros(1,nb-1);
		cU = -hankel(am); cV = hankel(bm,bp).';
	else
		bm = bp(2:end);
		am = am(2:end); ap = [am(end),zeros(1,nb-2)];
		cU = -hankel(am,ap);
		cV = hankel(bm).';
	end
end

% clean and resize R
[ cU, cV ] = compress_qr( cU, cV );




