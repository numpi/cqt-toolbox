function [cm,cp,cU,cV,cW,cZ]=fsi_tmult2(am, ap, bm, bp, m, p, n)
% function [cm,cp,cU,cV,cW,cZ]=fsi_tmult2(am, ap, bm, bp,n)
% Compute the product C=AB between an m x p  and an p x n Toeplitz matrices
% A and B defined by their first column and first row am,ap
% and bm,bp, respectively
% the matrix product C is written as toep(cm,cp) + cU*cV + 
% cW * cZ'(end:-1:1,end:-1:1)

lancz_param = 200; % size at which the lanczos method for 
                   % compressing is triggered
lancz_debug = 0;   % check the residual of the compression

if isempty(am) || isempty(bm)
    cm = [];
    cp = [];
    cU = [];
    cV = [];
    cW = [];
    cZ = [];
    return;
end

% Make sure all the inputs are row vectors
if ~isvector(am) || ~isvector(ap) || ~isvector(bm) || ~isvector(bp)
    error('The symbols must be specified as vectors');
end

am = reshape(am, 1, length(am));
ap = reshape(ap, 1, length(ap));
bp = reshape(bp, 1, length(bp));
bm = reshape(bm, 1, length(bm));

nam = length(am); nap = length(ap);
nbm = length(bm); nbp = length(bp);

% resize the input vectors according to their numerical length
[nam_eps,am] = epslength(am);
[nap_eps,ap] = epslength(ap);
[nbm_eps,bm] = epslength(bm);
[nbp_eps,bp] = epslength(bp);

% Convolution c(z) of a(z) e b(z), con calcolo del centro
a=[am(end:-1:1),ap(2:end)]; centroa=length(am);
b=[bm(end:-1:1),bp(2:end)]; centrob=length(bm);
c = conv_fft(a,b); cc = centroa+centrob-1;
cm = c(cc:-1:1);cp = c(cc:end);

% clean cm and cp
cm = truncate(cm,m);
cp = truncate(cp,n);
nam=length(am); nap=length(ap); nbm=length(bm); nbp=length(bp);


% Compute upper-left correction
h=min(nam,nbp)-1;

if h <= 0
    cU=[]; cV=[];
elseif h<=lancz_param % If at least one of the hankel 
                      % matrices is small then we compress with QR
    
    [ cU,cV ] = compress_qr(...
        hankel(-am(2:end),[-am(end),zeros(1,h-1)]),...
        hankel(bp(2:h+1),[bp(h+1:end),zeros(1,h-1)]).' );
elseif h>lancz_param   % If both the hankel matrices are big then 
                       % we compress with Lanczos
    
    cU=-am(2:end);
    cV=bp(2:end);
    if lancz_debug
        hh = min(length(cU),length(cV));
        Hu = hankel(cU); Hu = Hu(:,1:hh);
        Hv = hankel(cV); Hv = Hv(1:hh,:);
        H = Hu * Hv;
        clear Hu, Hv;
    end
    
    [ cU,cV ] =hankel_compress(cU, cV, cqtoption('compression'));
    
    if lancz_debug
        disp('*****residue upper-left corner*****')
        norm(cU*cV'-H)/norm(H)
        rank(H)
    end
    cV=conj(cV);
end


% Compute lower-right correction
if m >= p   % left factor
    if 1+m-p <= nam
        new_am = am(1+m-p:end);
    else
        new_am = 0;
    end
    new_ap = [zeros(1,max(0,1+m-p-nam)),am(min(1+m-p,nam):-1:2),ap];
else
    if 1+p-m <= nap
        new_ap = ap(1+p-m:end);
    else
        new_ap = 0;
    end
    new_am = [zeros(1,max(0,1+p-m-nap)),ap(min(1+p-m,nap):-1:2),am];
end
if p >= n   % right factor
    if 1+p-n <= nbm
        new_bm = bm(1+p-n:end);
    else
        new_bm = 0;
    end
    new_bp = [zeros(1,max(0,1+p-n-nbm)),bm(min(1+p-n,nbm):-1:2),bp];
else
    if 1+n-p <= nbp
        new_bp = bp(1+n-p:end);
    else
        new_bp = 0;
    end
    new_bm = [zeros(1,max(0,1+n-p-nbp)),bp(min(1+n-p,nbp):-1:2),bm];
end
am = new_am; ap = new_ap; bm = new_bm; bp = new_bp;
nam = length(am); nap = length(ap);
nbm = length(bm); nbp = length(bp);
h = min(nap,nbm)-1;
if h == 0
    cW=[]; cZ=[];
elseif h<=lancz_param && h>0 % If at least one of the hankel 
                             % matrices is small then we compress with QR
    [ cW,cZ ] = compress_qr(hankel(-ap(2:end),[-ap(end),zeros(1,h-1)]),...
        hankel(bm(2:h+1),[bm(h+1:end),zeros(1,h-1)]).' );
elseif h>lancz_param % If both the hankel matrices are big then 
                     % we compress with Lanczos
    cW=-ap(2:end);
    cZ=bm(2:end);
    if lancz_debug
        hh = min(length(cW),length(cZ));
        Hu = hankel(cW); Hu = Hu(:,1:hh);
        Hv = hankel(cZ); Hv = Hv(1:hh,:);
        H = Hu * Hv;
        clear Hu, Hv;
    end
    [ cW,cZ ] = hankel_compress(cW, cZ, cqtoption('compression')); 
    
    if lancz_debug
        disp('*****residue lower-right corner*****')
        norm(cW*cZ'-H)/norm(H)
        rank(H)
    end
    cZ=conj(cZ);
end



