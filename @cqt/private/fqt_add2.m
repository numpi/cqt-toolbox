%   Compute the sum of two finite quasi-Toeplitz matrices
function [ cm, cp, cu, cv, cw, cz ] = fqt_add2(am, ap, au, ...
    av,aw,az, bm, bp, bu, bv,bw,bz)
% sum the Toeplitz part
nam = length(am); nap = length(ap);
nbm = length(bm); nbp = length(bp);
cm = zeros(1,max(nam,nbm));
cp = zeros(1,max(nap,nbp));
cm(1:nam) = am; cm(1:nbm) = cm(1:nbm)+bm;
cp(1:nap) = ap; cp(1:nbp) = cp(1:nbp)+bp;

% sum the corrections
nau=size(au,1); nbu=size(bu,1); nu=max(nau,nbu);
nav=size(av,1); nbv=size(bv,1); nv=max(nav,nbv);
mau=size(au,2); mbu=size(bu,2);
mav=size(av,2); mbv=size(bv,2);

naw=size(aw,1); nbw=size(bw,1); nw=max(naw,nbw);
naz=size(az,1); nbz=size(bz,1); nz=max(naz,nbz);
maw=size(aw,2); mbw=size(bw,2);
maz=size(az,2); mbz=size(bz,2);
% compute U
u=zeros(nu,mau+mbu);
u(1:nau,1:mau)=au(1:nau,1:mau);
u(1:nbu,mau+1:mau+mbu)=bu(1:nbu, 1:mbu);
% compute V
v=zeros(nv,mav+mbv);
v(1:nav,1:mav)=av(1:nav,1:mav);
v(1:nbv,mav+1:mav+mbv)=bv(1:nbv, 1:mbv);
% compute W
w=zeros(nw,maw+mbw);
w(1:naw,1:maw)=aw(1:naw,1:maw);
w(1:nbw,maw+1:maw+mbw)=bw(1:nbw, 1:mbw);
% compute Z
z=zeros(nz,maz+mbz);
z(1:naz,1:maz)=az(1:naz,1:maz);
z(1:nbz,maz+1:maz+mbz)=bz(1:nbz, 1:mbz);
% compress and clean
nrm = fqt_norm(cm, cp, u, v, w, z);
[cu,cv] = compress_qr(u,v,nrm);
[cw,cz] = compress_qr(w,z,nrm);
% cm=cln(cm); cp=cln(cp);
[cm, cp] = symbol_clean(cm, cp, nrm);
end

