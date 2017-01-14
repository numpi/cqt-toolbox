function [cm, cp, Fc, Gc]=qt_inv(am, ap, F, G)
% compute the inverse C=A^{-1} of a semi-infinite QT-matrix A
% A = T + F*G',  T=toep(am,ap)
% T is defined by the first column am and its first row ap.
% C is a QT-matrix in the form C = toep(cm,cp) + Fc * Gc'
% The algorithm relies on the UL (spectral) factorization T = U*L
% and on the inversion formula
% C = Li*Ui - Li*((Ui*F)/Y)*(G'*Li)*Ui, Li=inv(L), Ui=inv(U), T=UL
% where Y=I+G'*Li*Ui*F
% algorithm:
% 1- T=UL
% 2- Ui=inv(U), Li=inv(L); [cm, cp, Fh, Gh]=Li*Ui
% 3- G1=Li'*G;  F1=Ui*F
% 4- Y=I+G1'*F1; F2=F1/Y; F3=Li*F2;  G2=Ui'*G1
% 5- Fc=[Fh, F3]; Gc=[Gh, G2]
% June 2, By Dario A. Bini

%1- Compute the spectral factorization T = UL

% The triangular case can be handled without calling spectral()
if length(am) == 1 || length(ap) == 1
	vm = am / am(1);
	vp = ap;
else
	[vm,vp] = spectral(am, ap);
end

%2- Invert the triangular Toeplitz matrices L and U
vm2 = vm(1); vp1 = vp(1);

[linvm, linvp] = reciprocal(vm,vm2);
[uinvm, uinvp] = reciprocal(vp1,vp);

% Compute Li*Ui
[cm, cp, Fh, Gh] = si_tmult(linvm, linvp, uinvm, uinvp);

%3-  compute G1 and F1
[nf,kf] = size(F); [ng, kg]=size(G);
[nfh,kfh] = size(Fh); [ngh,kgh]=size(Gh);
niu = length(uinvp); nil = length(linvm);
nf1 = nf; nf2=nf1; nf3 = nf2+nil-1;
ng1 = ng; ng2 = ng+niu-1;
F1 = toepmult_fft(uinvm, uinvp, nf1, nf, F);
G1 = toepmult_fft(linvp, linvm, ng1, ng, G);
%4- Compute Y, F2, G2 and F3
nfg=min(nf1,ng1);
Y = eye(kf)+G1(1:nfg,:).'*F1(1:nfg,:);
F2 = F1/Y; F3 = toepmult_fft(linvm, linvp, nf3, nf2, F2);
G2 = toepmult_fft(uinvp,uinvm,ng2,ng1 ,G1);
nf = max(nfh,nf3); ng = max(ngh,ng2);
Fc = zeros(nf,kf+kfh); Gc = zeros(ng,kf+kfh);
Fc(1:nfh,1:kfh) = Fh; Fc(1:nf3,kfh+1:kf+kfh) = -F3;
Gc(1:ngh,1:kfh) = Gh; Gc(1:ng2,kfh+1:kf+kfh) = G2;
[Fc,Gc] = compress_qr(Fc,Gc);


