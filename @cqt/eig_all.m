function [x, xcont, res, it, ei] = eig_all(AA, varargin)
% function [x1, xcont, res, it, ei] = eig_all(A, varargin)
% Compute the eigenvalue of the QT matrix A by functional iteration
% starting from the eigenvalues of a finite section of A of suitable size 
% In Input:
% A : the QT-matrix A=cqt(am,ap,E); 
% Optional parameters
% 'algo', a,  where a = 1,2,3,4, according to the algorithm used:
%      1 Newton's iteration applied to det(WV)=0 (Vandermonde version)
%      2 Newton's iteration applied to det(W[I;G;G^2;...)=0 (Frobenius version)
%      3 Newton's iteration applied to det(HV-x V)=0 (Vandermonde version)
%      4 Newton's iteration applied to det(H(G)-x I)=0 (Frobenius version)
%      default value: 1
% 'maxit', m, where m is the maximum number of iterations, default m=20
% 'epsilon', ep, where ep is the relative precision for the halt criterion,
%      default ep=1.e3*eps
% 'fact', f, the value f determines the size N=f*(max input size) of the matrix, default f=3
%      A_N whose eigenvalues are the initial points of the iterations, default 3
% 'verbose', ver, if ver=true, some information is printed at run time,
%      default ver=false
% 'plotfig', plf, if plf=true, some figures are plotted, default false
% 'advpx', ad, if ad=true, the multiprecision toolbox advanpix is used, 
%      default ad=false
% 'digits', dig, where dig is the number of decimal digits precision used in 
%      advanpix, default  dig=34
%
% In Output
% x: vector with the isolated eigenvalues  (p=q)
% xcont: eigenvalues in a continuous set
% res: vector with the residual errors of x
% it: vector with the number of iterations
% ei: vector with the eigenvalues of the finite section A_N


% By D.A. Bini, L. Robol, March,4 , 2022

p = inputParser;
addOptional(p, 'algo', 1);
addOptional(p, 'maxit', 20);
addOptional(p, 'epsilon', []); % Chosen based on advpx
addOptional(p, 'fact', 3);
addOptional(p, 'verbose', false);
addOptional(p, 'plotfig', false);
% addOptional(p, 'residual', false);
addOptional(p, 'advpx', false);
addOptional(p, 'digits', 34);

parse(p, varargin{:});

algo = p.Results.algo;
maxit = p.Results.maxit;
epsi = p.Results.epsilon;
fact = p.Results.fact;
verbose = p.Results.verbose;
plotfig = p.Results.plotfig;
% residual = p.Results.residual;
advpx = p.Results.advpx;
digits = p.Results.digits;

if isempty(epsi)
    if advpx
        epsi = mp(10)^(-digits+3);
    else
        epsi = 1e3 * eps;
    end
end
[am, ap] = symbol(AA);
E = correction(AA);

% adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end

% check for trailing zeros in am, ap
  nm = length(am); np = length(ap);
  k= 0;
  for i=nm:-1:1
    if am(i)~=0
      break
    else
      k = k+1;
    end
  end
  nm = nm-k; m = nm-1;  
  k= 0;
  for i=np:-1:1
    if ap(i)~=0
      break
    else
      k = k+1;
    end
  end
  np = np-k; n = np-1;  
  am = am(1:nm); ap = ap(1:np); AA = cqt(am,ap,E);


% compute initial approximations
  n = max([length(am) + length(ap),size(E,1),size(E,2)]); %%% 4/3/22
  n = n*fact;
  if verbose
      fprintf('Computing the eigenvalues of the truncated matrix of size %d\n',n);
  end
  A=single(toepl(am,ap,E,n));
  tic;  ei = eig(A);  teig = toc;
  ei = double(ei);
  if advpx
     mp.Digits(digits);
     ei = mp(ei);
     am = mp(am); ap = mp(ap);
  end

  tic;
  c1 = 1; c2 = 1; 
  x = []; xcont = []; res=[]; it=[];
  if advpx
     x = mp(x); xcont = mp(xcont); res = mp(res);
  end
  if am(1)~=ap(1)
      fprintf('WARNING: am(1) and ap(1) have inconsistent values\n')
      return
  end

  for j=1:n
     w = wind(am,ap,ei(j),advpx);
     if w>=0  %%%%%%%%    
         [z,~,~,info]=eig_single(AA, ei(j), 'algo', algo, 'maxit', maxit, ...
                'epsilon', epsi, 'verbose', verbose, ...
                'residual', false, 'advpx', advpx, 'digits', digits);
         if info.flg == 1 || info.flg == 3         
            x(c1) = z;
            it(c1) = info.it;
            c1=c1+1; 
         elseif info.flg == 2
            xcont(c2) = z; c2 = c2+1;
         end
     end
  end
  tqt = toc;
  
  if  ~isempty(x)      % refine the approximations and compute the residual errors   
    maxit = 1; epsix=1.e-20;
    if verbose
      fprintf('\n eig_all: refining the approximations....\n');
    end
    for j=1:length(x)
        [x(j),~,~,info]=eig_single(AA, x(j), 'algo', algo, 'maxit', maxit, ...
              'epsilon', epsi*100, 'verbose', verbose, ...
              'residual', true, 'advpx', advpx, 'digits', digits);
        res(j) = info.res;
    end
% remove copies of the same eigenvalue from x1
    [x ,res, it] = clean(x,res,it,epsi);
  end

% plot figures
  if plotfig
     figure; 
     plot(ei, 'ro');  hold on;
     plot(x+1.e-200*1i,'b.', 'markersize',10);
     ax = gca; ax.FontSize = 16;
     plot(xcont+1.e-20*1i, '.c');
     range(AA);
     drawnow;
     hold off;
  end
  if verbose
     fprintf('Time required by eig = %e time required by eig_all =  %e\n',teig,tqt); 
  end
end



function [y,rs,iter] = clean(x,res,it,epsi)
% function y = clean(x,epsi)
% this function removes copies of the same component in x
% and keeps the values having the minimum residual

  n = length(x);
  [xs,xi] = sort(res,'ascend');
  res = xs;
  x = x(xi); it = it(xi);
  ax = abs(x);
  y(1)=x(1); 
  rs(1) = res(1); 
  iter(1) = it(1);  c = 1;
  for j=2:n
     flg = true;
     for k=1:c
        if abs(x(j) - y(k))<epsi*ax(j)*10
           flg = false;
           break
        end
     end
     if flg
        c=c+1;
        y(c) = x(j); iter(c) = it(j);         
        rs(c) = res(j); 
     end
  end
% sort
  [~,idx] = sort(real(y));
  y = y(idx); iter =iter(idx); 
  rs = rs(idx);
end

function a = toepl(am,ap,E,n)
% function a = toepl(am,ap,E,n)
% The nxn matrix a = CQT(am,ap,E) is computed
nm = length(am); np = length(ap);
%n = max(n,(nm+np)*2);  
a = am(1)*eye(n);
for j=2:nm
  a = a+am(j)*diag(ones(n-j+1,1),-j+1);
end
for j=2:np
  a = a+ap(j)*diag(ones(n-j+1,1),j-1);
end
[k1,k2]=size(E);
a(1:k1,1:k2) = a(1:k1,1:k2)+E;
end
