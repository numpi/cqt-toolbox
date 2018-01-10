% Compute C=A^{-1} of a semi-inUinite QT-matrix A
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
function [cm, cp, Uc, Vc,Wc,Zc]=fqt_inv2(am, ap, aU, aV, aW,aZ,n)
%1- Compute the spectral factorization T = UL
if am(1) == 0
	error('fqt_inv2:: Can not invert a symbol with a0 = 0');
end
% The triangular cases can be handled without calling spectral()
if length(am) == 1 || length(ap) == 1
	vm = am / am(1);
	vp = ap;
else
	[vm, vp] = spectral(am,ap);
end

% Compute the lower right corner correction
[~,~,~,~,hW,hZ] = fsi_tmult2(vp(1), vp, vm, vm(1),n, n, n);
m1=size(aW,1); m2=size(hW,1);
if (m1 < m2)
	aW = [aW;zeros(m2-m1,size(aW,2))];
elseif(m2 < m1)
	hW = [hW;zeros(m1-m2,size(hW,2))];
end
aW = [aW,-hW];
m1 = size(aZ,1); m2 = size(hZ,1);
if (m1 < m2)
	aZ = [aZ;zeros(m2-m1,size(aZ,2))];
elseif(m2 < m1)
	hZ = [hZ;zeros(m1-m2,size(hZ,2))];
end
aZ = [aZ,hZ];
[aW,aZ] = compress_qr(aW,aZ);

%2- Invert the triangular Toeplitz matrices L and U
[linvm, linvp] = reciprocal(vm,vm(1), n, n);
[uinvm, uinvp] = reciprocal(vp(1),vp, n, n);
linvm = truncate(linvm,n);  uinvp = truncate(uinvp,n);

% Compute Li*Ui
[cm, cp, cU, cV] = fsi_tmult2(linvm, linvp, uinvm, uinvp, n, n, n);

% 3-  compute G1 and F1

[nU,kU] = size(aU); [nV] = size(aV,1);
[nW,kW] = size(aW); [nZ] = size(aZ,1);
niu = length(uinvp); nil = length(linvm);
[nUh,kUh] = size(cU); [nVh]=size(cV,1);
tdim1=min(nU+nil-1,n); tdim2=min(n,niu+nW-1);
tdim3=min(n,nil+nZ-1); tdim4=min(nV+niu-1,n);
aU1 = toepmult_fft(uinvm, uinvp, nU, nU, aU);
aV1 = toepmult_fft(linvp, linvm, nV, nV, aV);
aW1 = toepmult_fft(uinvp, uinvm, tdim2, nW, aW);
aZ1 = toepmult_fft(linvm, linvp, tdim3, nZ, aZ);

%4- Compute Y, F2, G2 and F3
nU1 = size(aU1,1); nV1=size(aV1,1);
nW1 = size(aW1,1); nZ1=size(aZ1,1);
Y = eye(kU + kW);
nmin = min(nU1,nV1);
Y(1:kU,1:kU) = Y(1:kU,1:kU) + aV1(1:nmin,:).' * aU1(1:nmin,:);
nmin = min(size(aW1,1),size(aZ1,1));
temp = aZ1(1:nmin,:).' * aW1(1:nmin,:);
Y(end-kW+1:end,end-kW+1:end)=Y(end-kW+1:end,end-kW+1:end)+temp(end:-1:1,end:-1:1);


if(nV1+nW1<=n && nU1+nZ1<=n)
	aU2 = aU1/Y(1:kU,1:kU);
	Y = Y(end:-1:1,end:-1:1);
	aW2 = aW1/Y(1:kW,1:kW);
else
	if(nV1+nW1>n)
		Y(1:kU,kU+1:end) = Y(1:kU,kU+1:end)+aV1(n-nW1+1:nV1,1:kU).'*aW1(end:-1:end-(nV1+nW1-n)+1,end:-1:1);
	end
	if(nU1+nZ1>n)
		Y(kU+1:end,1:kU) = Y(kU+1:end,1:kU)+aZ1(end:-1:end-(nU1+nZ1-n)+1,end:-1:1).'*aU1(n-nZ1+1:nU1,1:kU);
	end
	Y = inv(Y);
	Y = [[aU1;zeros(n-nU1,kU)],[zeros(n-nW1,kW);aW1(end:-1:1,end:-1:1)]]*Y;
	aU2 = Y(:,1:kU);
	aW2 = Y(:,kU+1:end);
	aW2 = aW2(end:-1:1,end:-1:1);
	
end
nU = size(aU2,1);
nW = size(aW2,1);
tdim1 = min(nU+nil-1,n); tdim2 = min(n,niu+nW-1);
aU3 = toepmult_fft(linvm, linvp, tdim1, nU, aU2);
aW3 = toepmult_fft(linvp, linvm, nW, nW, aW2);
aV2 = toepmult_fft(uinvp,uinvm,tdim4,nV ,aV1);
aZ2 = toepmult_fft(uinvm,uinvp, tdim3,tdim3 ,aZ1);

% 5- Fc=[Fh, F3]; Gc=[Gh, G2]
nUmax = max(nUh,tdim1); nVmax = max(nVh,tdim4);
Uc = zeros(nUmax,kU+kUh); Vc = zeros(nVmax,kU+kUh);
Uc(1:nUh,1:kUh) = cU;
Uc(1:tdim1,kUh+1:kU+kUh) = -aU3;
Vc(1:nVh,1:kUh) = cV; Vc(1:tdim4,kUh+1:kU+kUh) = aV2;

Wc = -aW3;
Zc = aZ2;

nrm = fqt_norm(cm, cp, Uc, Vc, Wc, Zc);
[cm, cp] = symbol_clean(cm, cp, nrm);
[Uc,Vc] = compress_qr(Uc,Vc, nrm);
[Wc,Zc] = compress_qr(Wc,Zc, nrm);
