function [ cm, cp, cu, cv ] = qt_add(am, ap, au, av, bm, bp, bu, bv)
%QT_ADD Sum two CQT matrices.
%
%   Sum of two quasi-Toeplitz matrices T(am,ap)+au*av', T(bm,bp)+bu*bv'
%   The sum is T(cm,cp)+cu*cv'
% June 2, By Dario A. Bini

% sum the Toeplitz part
nam = length(am); nap = length(ap);
nbm = length(bm); nbp = length(bp);
cm = zeros(1,max(nam,nbm));
cp = zeros(1,max(nap,nbp));
cm(1:nam) = am; cm(1:nbm) = cm(1:nbm)+bm;
cp(1:nap) = ap; cp(1:nbp) = cp(1:nbp)+bp;

% sum the corrections
%  detect the sizes
nau = size(au, 1); nbu = size(bu, 1); nu = max(nau, nbu);
nav = size(av, 1); nbv = size(bv, 1); nv = max(nav, nbv);
mau = size(au, 2); mbu = size(bu, 2);
mav = size(av, 2); mbv = size(bv, 2);
% compute U
u = zeros(nu,mau+mbu);
u(1:nau, 1:mau) = au(1:nau, 1:mau);
u(1:nbu, mau+1:mau+mbu) = bu(1:nbu, 1:mbu);
% compute V
v = zeros(nv, mav+mbv);
v(1:nav, 1:mav) = av(1:nav, 1:mav);
v(1:nbv, mav+1:mav+mbv) = bv(1:nbv, 1:mbv);
% compress and clean
nrm = qt_norm(cm, cp, u, v);
[cu, cv] = compress_qr(u, v, nrm);
[cm, cp] = symbol_clean(cm, cp, nrm);


