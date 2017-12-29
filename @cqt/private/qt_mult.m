function [cm,cp,cU,cV]=qt_mult(am, ap, aU,aV, bm, bp, bU,bV)
%QT_MULT Compute the product of two CQT matrices.
%
% Computes the product C=AB of two semi infinite
% quasi toeplitz matrices A=toep(am,ap)+Ua*Va', B=toep(bm,bp)+Ub*Vb'
% am, ap, and bm, bp are the first column and the first row
% of the Toeplitz part of A and B, respectively
% the vectors cm and cp are first column and first row of
% C=toep(cm,cp) + Uc*Vc'
% [Uc, Vc] = toep(am,ap)*Ub*Vb' + Ua*Va'*toep(bm,bp) + Ua(Va'Ub)*Vb + R
% where  R: toep(a)*toep(b)=toep(c) + R, R=Ur*Vr'
% June 2, By Dario A. Bini

% Start the computation
% compute the products
% 0  [cm,cp,rU,rV]= si_tmult(am,ap,bm,bp);
% 1  TRb = toepmul(am,ap,rc1,rb1)*Rb;
% 2  RaT = Ra*toep(bm,bp,ra2,rc2);
% 3  RaRb = Ra(:,1:min(ra2,rb1))*Rb(1:min(ra2,rb1),:);

% 0 toeplitz product
[cm, cp, rU, rV] = si_tmult(am, ap, bm, bp);
ra1 = size(aU, 1); ra2 = size(aV, 1); ka = size(aU, 2);
rb1 = size(bU, 1);rb2 = size(bV, 1); kb = size(bU, 2);
r1 = size(rU, 1); r2 = size(rV, 1); kr = size(rU, 2);
c1 = length(cm); c2 = length(cp);
a1 = length(am); a2 = length(ap);
b1 = length(bm); b2 = length(bp);
tc1 = rb1 + a1;    %dim1 of toep(a)*Rb
tc2 = ra2 + b2;    %dim2 of Ra*toep(b)
% compute max dim of Rc
rc1 = max([r1, tc1, ra1]); %dim1 of R in output
rc2 = max([rb2, r2, tc2]);
% 1
TRbU = toepmult_fft(am,ap,rc1,rb1,bU) ; TRbV = bV;   % 1
% 2
RaTU = aU; RaTV = toepmult_fft(bp,bm,rc2,ra2,aV);
sav1 = size(aV,1); sbu1 = size(bU,1);
mx = max(sav1,sbu1);
kav = size(aV,2); kbu = size(bU,2);
VV = zeros(mx,kav); UU = zeros(mx,kbu);
VV(1:sav1,:) = aV;
UU(1:sbu1,:) = bU;
%  3
RaRbU = aU*(VV.'*UU); RaRbV = bV;
% sum them up:  cU = [rU, TRbU, RaTU, RaRbU]; cV=[rV, TRbV, RaTV,RaRbV];
nru = size(rU,1); ntrbu = size(TRbU,1);
nratu = size(RaTU,1); nrarbu = size(RaRbU,1);
nrv = size(rV,1); ntrbv = size(TRbV,1);
nratv = size(RaTV,1); nrarbv = size(RaRbV,1);
kru = size(rU,2); ktrbu = size(TRbU,2);
kratu = size(RaTU,2); krarbu = size(RaRbU,2);
rc1 = max([nru,ntrbu,nratu,nrarbu]);
rc2 = max([nrv,ntrbv,nratv,nrarbv]);
k = kru+ktrbu+kratu+krarbu;
cU = zeros(rc1,k);
cV = zeros(rc2,k);
cU(1:nru,1:kru) = rU;
cV(1:nrv,1:kru) = rV;
cU(1:ntrbu,kru+1:kru+ktrbu) = TRbU;
cV(1:ntrbv,kru+1:kru+ktrbu) = TRbV;
cU(1:nratu,kru+ktrbu+1:kru+ktrbu+kratu) = RaTU;
cV(1:nratv,kru+ktrbu+1:kru+ktrbu+kratu) = RaTV;
cU(1:nrarbu,kru+ktrbu+kratu+1:kru+ktrbu+kratu+krarbu) = RaRbU;
cV(1:nrarbv,kru+ktrbu+kratu+1:kru+ktrbu+kratu+krarbu) = RaRbV;
% compress and clean
nrm = qt_norm(cm, cp, cU, cV);
[cU, cV] = compress_qr(cU, cV, nrm);
[cm, cp] = symbol_clean(cm, cp, nrm);



