function [cm, cp, Fc, Gc] = qt_mldivide(am, ap, F, G, bm, bp, W, Z)
% compute the solution of the linear system A x = B where A is a cqt-matrix

% 0 - [L,U] = spectral factorization of toep(am,ap)
% 1 - U^-1 W
% 2 - L^-1 U^-1 W
% 3- G L^-1 U^-1 W
mW = size(W,1);
mF = size(F,1);
switch cqtoption('inversion')
    case 'cr'
        spectral = @spectral_cr;
        reciprocal = @reciprocal_cr;
    case 'fft'
        spectral = @spectral_fft;
        reciprocal = @reciprocal_fft;
end

%1- Compute the spectral factorization T = UL
  [vm,vp] = spectral(am,ap);

%2- Invert the triangular Toeplitz matrices L and U

  

  vm2 = [vm(1)]; vp1 = [vp(1)];  
  [linvm, linvp] = reciprocal(vm,vm2);
  [uinvm, uinvp] = reciprocal(vp1,vp);

  mG1 = min(length(linvm) + mW, size(G, 1));
  mG2 = min(length(linvm) + mF, size(G, 1));

  [cm,cp,cU,cV] = qt_mult(uinvm, uinvp, 0, 0, bm, bp, 0, 0);
  [cm,cp,cU,cV] = qt_mult(linvm, linvp, 0, 0, cm,cp,cU,cV);
  iUW = toepmult_fft(uinvm, uinvp, mW, mW, W);
  iLiUW = toepmult_fft(linvm, linvp, length(linvm) + mW, mW, iUW);
  iUF = toepmult_fft(uinvm, uinvp, mF, mF, F);
  iLiUF = toepmult_fft(linvm, linvp, length(linvm) + mF, mF, iUF);
% 3- 
 GiLiUW = G(1:mG1 ,:).'*iLiUW(1:mG1 ,:);
 GiLiUF = G(1:mG2 ,:).'*iLiUF(1:mG2 ,:);

% Computation of the central block of Woodbury identity
S = eye(size(F, 2)) + GiLiUF;

% 4-
mG3 = min(size(cU,1), size(G, 1));
[tiLiUW,tiLiUF] = align_mat( iLiUW,iLiUF,'rows');
tFc = tiLiUW - tiLiUF * (S \  GiLiUW);
tFFc = iLiUF * (S \  G(1:mG3,:).' * cU(1:mG3,:));

[FFc1, FFc2] = align_mat(cU, tFFc, 'both');
FFc = FFc1 - FFc2;

[Fc1, Fc2] = align_mat(tFc, FFc, 'rows');
Fc = [Fc1, Fc2];

[Gc1, Gc2] = align_mat(Z, cV, 'rows');
Gc = [Gc1, Gc2];

Gc3 = toepmult_fft(cp, cm, size(G,1) + length(cp), size(G,1), G);
Fc3 = -iLiUF / S;

[Gc,Gc3] = align_mat(Gc, Gc3, 'rows');
Gc = [ Gc, Gc3 ];
[Fc, Fc3] = align_mat(Fc, Fc3, 'rows');
Fc = [Fc Fc3];

[Fc,Gc] = compress_qr(Fc,Gc);
