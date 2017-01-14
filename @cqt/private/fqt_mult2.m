function [cm,cp,cU,cV,cW,cZ]=fqt_mult2(am, ap, aU, aV, aW, aZ, bm, bp, bU, bV, bW, bZ, m, p, n)
% Computes the product C=AB between an m x p and an p x n   
% quasi toeplitz matrices A=toep(am,ap)+Ua*Va' + Wa *Za'(end:-1:1,end:-1:1), B=toep(bm,bp)+Ub*Vb' + Wb *Zb'(end:-1:1,end:-1:1)
% am, ap, and bm, bp are the first column and the first row
% of the Toeplitz part of A and B, respectively
% the vectors cm and cp are first column and first row of 
% C=toep(cm,cp) + Uc*Vc'
% [Uc, Vc] = toep(am,ap)*Ub*Vb' + Ua*Va'*toep(bm,bp) + Ua(Va'Ub)*Vb + R
% where  R: toep(a)*toep(b)=toep(c) + R, R=Ur*Vr'
%
% compute the products
% 0  [cm,cp,rU,rV,sW,sZ]= fsi_tmult2(am,ap,bm,bp);
% 1  TRb = toepmul(am,ap,rc1,rb1)*Rb; and T * Sb
% 2  RaT = Ra*toep(bm,bp,ra2,rc2); and Sa * T
% 3  RaRb = Ra(:,1:min(ra2,rb1))*Rb(1:min(ra2,rb1),:); and Sa *Sb
% 4  if the corrections are too big compute also RaSb e SaRb

% toeplitz product
  [cm,cp,rU,rV,sW,sZ] = fsi_tmult2(am, ap, bm, bp, m, p, n);   % rU,rV % 0
%  [cm,cp,rU,rV,sW,sZ] = x(am,ap,bm,bp,n);   % rU,rV % 0
  ra1 = size(aU,1); ra2 = size(aV,1); kar = size(aU,2);
  rb1 = size(bU,1);rb2 = size(bV,1); kbr = size(bU,2);
  r1 = size(rU,1); r2 = size(rV,1); kr = size(rU,2);

  sa1=size(aW,1); sa2=size(aZ,1); kas=size(aW,2);
  sb1=size(bW,1);sb2=size(bZ,1); kbs=size(bW,2);
  s1=size(sW,1); s2=size(sZ,1); ks=size(sW,2);

  c1 = length(cm); c2 = length(cp);
  a1 = length(am); a2 = length(ap);
  b1 = length(bm); b2 = length(bp);

% Symbols corresponding to the bottom right correction
if m >= p   % left factor
	if 1+m-p <= a1
		down_am = am(1+m-p:end);
	else
		down_am = 0;
	end
	down_ap = [zeros(1,max(0,1+m-p-a1)),am(min(1+m-p,a1):-1:2),ap];
else
	if 1+p-m <= a2
		down_ap = ap(1+p-m:end);
	else
		down_ap = 0;
	end
	down_am = [zeros(1,max(0,1+p-m-a2)),ap(min(1+p-m,a2):-1:2),am];
end
if p >= n   % right factor
	if 1+p-n <= b1
		down_bm = bm(1+p-n:end);
	else
		down_bm = 0;
	end
	down_bp = [zeros(1,max(0,1+p-n-b1)),bm(min(1+p-n,b1):-1:2),bp];
else
	if 1+n-p <= b2
		down_bp = bp(1+n-p:end);
	else
		down_bp = 0;
	end
	down_bm = [zeros(1,max(0,1+n-p-b2)),bp(min(1+n-p,b2):-1:2),bm];
end
da1 = length(down_am); da2 = length(down_ap);
db1 = length(down_bm); db2 = length(down_bp);

  tc1 = min(rb1 + a1,m);    %dim1 of toep(a)*Rb
  tc2 = min(ra2 + b2,n);    %dim2 of Ra*toep(b)

  ts1 = min(sb1 + da2,m);    %dim1 of toep(a)*Sb
  ts2 = min(sa2 + db1,n);    %dim2 of Sa*toep(b)



% compute max dim of R and S
  rc1 = max([r1,tc1,ra1]); %dim1 of R in output
  rc2 = max([rb2,r2,tc2]);

  sc1 = max([s1,ts1,sa1]); %dim1 of S in output
  sc2 = max([sb2,s2,ts2]);



  TRbU=toepmult_fft(am,ap,rc1,rb1,bU) ; TRbV = bV;   % 1
  TSbW=toepmult_fft(down_ap,down_am,sc1,sb1,bW) ; TSbZ = bZ;

  RaTU=aU; RaTV = toepmult_fft(bp,bm,rc2,ra2,aV);  % 2
  sav1 = size(aV,1); sbu1 = size(bU,1);
  mx = max(sav1,sbu1);
  kav = size(aV,2); kbu = size(bU,2);
  VV=zeros(mx,kav); UU = zeros(mx,kbu);
  VV(1:sav1,:) = aV; 
  UU(1:sbu1,:) = bU;

  SaTW = aW; SaTZ = toepmult_fft(down_bm,down_bp,sc2,sa2,aZ);  
  saz1 = size(aZ,1); sbw1 = size(bW,1);
  mx=max(saz1,sbw1);
  kaz = size(aZ,2); kbw = size(bW,2);
  ZZ = zeros(mx,kaz); WW = zeros(mx,kbw);
  ZZ(1:saz1,:) = aZ; 
  WW(1:sbw1,:) = bW;

% RaRbU=aU*(aV'*bU); RaRbV = bV;                   % 3
  RaRbU = aU*(VV.'*UU); RaRbV = bV;
  SaSbW = aW*(ZZ.'*WW); SaSbZ = bZ;

% Case of too big corrections ---> we need to consider cross products
 krasb = 0;
 ksarb = 0;
 up1=0; up2=0; down1=0; down2=0;
  if(ra2+sb1>p)
	krasb = ra2+sb1-p;
	WW = bW(end:-1:1,end:-1:1);
	rasbV = bZ(end:-1:1,end:-1:1);	
	rasbU = aU*(aV(end-krasb+1:end,:).'*WW(1:krasb,:));
	if(ra2>=sb2)  % We decide to locate the correction in the upper left or in the lower right corner
		up1 = 1;
	else
		down1 = 1;
	end
  end
  if(sa2+rb1>p)
	ksarb = sa2+rb1-p;
	sarbV = bV;
	WW = aW(end:-1:1,end:-1:1);
	ZZ = aZ(end:-1:1,end:-1:1);
	sarbU = WW*(ZZ(1:ksarb,:).'*bU(end-ksarb+1:end,:));
	if(rb2>=sa2)
		up2 = 1;
	else
		down2 = 1;
	end
  end	

 % sum them up % cU = [rU, TRbU, RaTU, RaRbU, rasbU,sarbU]; cV=[rV, TRbV, RaTV,RaRbV,rasbV, sarbV]
 % Compute the upper left correction
  nru = size(rU,1); ntrbu = size(TRbU,1); nratu = size(RaTU,1); nrarbu = size(RaRbU,1);
  nrv = size(rV,1); ntrbv = size(TRbV,1); nratv = size(RaTV,1); nrarbv = size(RaRbV,1);
  kru = size(rU,2); ktrbu = size(TRbU,2); kratu = size(RaTU,2); krarbu = size(RaRbU,2);
  rc1 = max([nru,ntrbu,nratu,nrarbu,up1*ra1,up2*m]);
  rc2 = max([nrv,ntrbv,nratv,nrarbv,up1*n,up2*rb2]);
  k = kru+ktrbu+kratu+krarbu+kbs*up1+ kbr*up2;
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
  
  if(up1==1)   % If the result of the cross products is located up
	cU(1:ra1,kru+ktrbu+kratu+krarbu+1:kru+ktrbu+kratu+krarbu+kbs) = rasbU;
	cV(end-sb2+1:end,kru+ktrbu+kratu+krarbu+1:kru+ktrbu+kratu+krarbu+kbs) = rasbV;
  end
  if(up2==1)
	cU(end-sa1+1:end,kru+ktrbu+kratu+krarbu+up1*kbs+1:kru+ktrbu+kratu+krarbu+up1*kbs+kbr) = sarbU;
	cV(1:rb2,kru+ktrbu+kratu+krarbu+up1*kbs+1:kru+ktrbu+kratu+krarbu+up1*kbs+kbr) = sarbV;
  end	
  % Compute the lower right correction
  nsw = size(sW,1); ntsbw = size(TSbW,1); nsatw = size(SaTW,1); nsasbw = size(SaSbW,1);
  nsz = size(sZ,1); ntsbz = size(TSbZ,1); nsatz = size(SaTZ,1); nsasbz = size(SaSbZ,1);
  ksw = size(sW,2); ktsbw = size(TSbW,2); ksatw = size(SaTW,2); ksasbw = size(SaSbW,2);
  sc1 = max([nsw,ntsbw,nsatw,nsasbw,down1*m,down2*sa1]);
  sc2 = max([nsz,ntsbz,nsatz,nsasbz,down1*sb2,down2*n]);
  k = ksw+ktsbw+ksatw+ksasbw+kbs*down1+ kbr*down2;
  cW = zeros(sc1,k);
  cZ = zeros(sc2,k);
  cW(1:nsw,1:ksw) = sW;
  cZ(1:nsz,1:ksw) = sZ;
  cW(1:ntsbw,ksw+1:ksw+ktsbw) = TSbW;
  cZ(1:ntsbz,ksw+1:ksw+ktsbw) = TSbZ;
  cW(1:nsatw,ksw+ktsbw+1:ksw+ktsbw+ksatw) = SaTW;
  cZ(1:nsatz,ksw+ktsbw+1:ksw+ktsbw+ksatw) = SaTZ;
  cW(1:nsasbw,ksw+ktsbw+ksatw+1:ksw+ktsbw+ksatw+ksasbw) = SaSbW;
  cZ(1:nsasbz,ksw+ktsbw+ksatw+1:ksw+ktsbw+ksatw+ksasbw) = SaSbZ;
  
  if(down1==1)    % If the result of the cross products is located down
	cW(end-ra1+1:end,ksw+ktsbw+ksatw+ksasbw+1:ksw+ktsbw+ksatw+ksasbw+kbs) = rasbU(end:-1:1,end:-1:1);
	cZ(1:sb2,ksw+ktsbw+ksatw+ksasbw+1:ksw+ktsbw+ksatw+ksasbw+kbs) = rasbV(end:-1:1,end:-1:1);
  end
  if(down2==1)
	cW(1:sa1,ksw+ktsbw+ksatw+ksasbw+down1*kbs+1:ksw+ktsbw+ksatw+ksasbw+down1*kbs+kbr) = sarbU(end:-1:1,end:-1:1);
	cZ(end-rb2+1:end,ksw+ktsbw+ksatw+ksasbw+down1*kbs+1:ksw+ktsbw+ksatw+ksasbw+down1*kbs+kbr) = sarbV(end:-1:1,end:-1:1);
  end
% compress and clean
[cU,cV] = compress_qr(cU,cV);
[cW,cZ] = compress_qr(cW,cZ);
  if (min(size(cU)==0) || min(size(cV)==0))
	cU = []; cV = [];
 
  end
  if (min(size(cW)==0) || min(size(cZ)==0))
	cW = []; cZ = [];
 
  end
  cm = cln(cm); cp = cln(cp);

