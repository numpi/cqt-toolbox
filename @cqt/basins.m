function [AA,v,iter] = basins(A,n,x0,x1,y0,y1,varargin) 
% function [A,v,iter]  = basins(A,n,x0,x1,y0,y1,varargin) 
% If plotfig is true then the function draws the basins of attraction of the 
% fixed point iteration given by algorithm algo, where algo=1,2,3,4,
% applied to the CQT matrix A in the range x0<re(x)<x1, y0<im(x)<y1
% 
% Optional parameters:
% 'algo', a,  where a = 1,2,3,4, according to the algorithm used:
%       1 Newton's iteration applied to det(WV)=0 (Vandermonde version)
%       2 Newton's iteration applied to det(W[I;G;G^2;...)=0 (Frobenius version)
%       3 Newton's iteration applied to det(HV-x V)=0 (Vandermonde version)
%       4 Newton's iteration applied to det(H(G)-x I)=0 (Frobenius version)
%       default value: 1
% 'maxit', m, where m is the maximum number of iterations, default m=20
% 'epsilon', ep, where ep is the relative precision for the halt criterion,
%        default ep=1.e-11
% 'verbose', ver, if ver=true, some information is printed at run time,
%        default ver=false
% 'advpx', ad, if ad=true, the multiprecision toolbox advanpix is used, 
%        default ad=false
% 'digits', dig, where dig is the number of decimal digits precision used in 
%        advanpix, default  dig=34
% 'plotfig', if true, the figure is plotted, default true

% On output: A is the nxn matrix of unint n: A(k,j) contains the integer s 
% such that the fixed point sequence obtained by starting from the point 
% z0 = x0+(k-1)(x1-x0)/(n-1) + i(y0+(j-1)(y1-y0)/(n-1)), k,j=1:n
% converges to the (s+1)-st eigenvalue. Moreover,
% A(k,j) = 256-1 if z0 is in a continuous set
%          255-1 if the sequence exits the component containing z0
%          254-1 if the sequence diverges or if exceeded the max number of iterations
%
% v:    vector with the isolated eigenvalues
% iter: overall number of iterations
% a picture of the basins is plotted with the following colors:
% x eig: color associated with the eig      (flag 1 or 3)
% 255 continuous set:       green           (flag 2)
% 254 out of the component: light gray      (flag 5)
% 253 p<q not an eig :      black           (flag 4)
% 253 Exceeded iter:        black           (flag 6)
% 253 Diverging seq.        black           (flag 7)

% By D.A. Bini, May 3, 2021

info = struct;

p = inputParser;%algo,maxit, epsilon,verbose,plotfig,advpx,digits)

addOptional(p, 'algo', 1);
addOptional(p, 'maxit', 20);
addOptional(p, 'epsilon', 1.e-11);
addOptional(p, 'verbose', false);
addOptional(p, 'plotfig', true);
addOptional(p, 'advpx', false);
addOptional(p, 'digits', 34);

parse(p, varargin{:});

algo = p.Results.algo;
maxit = p.Results.maxit;
epsi = p.Results.epsilon;
verbose = p.Results.verbose;
advpx = p.Results.advpx;
digits = p.Results.digits;
plotfig = p.Results.plotfig;

  [am, ap] = symbol(A);
  E = correction(A);
  AA = uint8(zeros(n));
  v = []; 
  
  % adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end

% if x0=x1 set default values for x0 and x1  
  if x0==x1
     a = [ap(end:-1:1);am(2:end)];
     mxval = sum(abs(a)) + norm(E,inf);
     x0 = -mxval*1.001; x1 = mxval; y0 = x0; y1 = x1;
  end
  
  % compute the eigenvalues 
  v = eig_all(A, 'algo',algo);
  nz = length(v);
  
  % choose the color palette
   MP = rand(256,3);
   MP(256,:) = [0 1 0];   %  continuous set
   MP(255,:) = [0.5,0.5,0.5];   % out of the component
   MP(254,:) = [0 0 0];         % non-converging sequence

% construct the basins  
  iter = 0;
  for k=1:n
    for j=1:n
       z0 = x0+(j-1)*(x1-x0)/(n-1) + 1i*(y0+(k-1)*(y1-y0)/(n-1));
       [z, ~,~,info]= eig_single(A,z0,'algo',algo,'maxit',maxit,'epsilon',epsi,'advpx',advpx,'digits',digits);
       flg = info.flg; it = info.it;
       iter = iter+it;
       if verbose && j==1 &&  mod(k,10)==1
          fprintf('[k,j,flg]=[%d,%d,%d] \n',k,j,flg);
       end 
       if flg==1  || flg==3      %  isolated eigenvalue
          for l=1:nz
             if abs(z-v(l))<1.e-8 
                AA(k,j) = min(l-1,nz-1); 
                break
             end
          end
       elseif flg==2    % p>q: continuous eigenvalue
          AA(k,j) = 256-1;
       elseif flg==5    % change of component
          AA(k,j) = 255-1;
       elseif flg>=4    % non-converging sequence 
          AA(k,j) = 254-1;
       end
    end
  end  
  if plotfig
     figure; image([x0,x1],[y0,y1],AA); %title(' Basins ');
     colormap(MP);
     hold on;
     range(A);
     hold off;
     ax = gca;
     ax.YDir = 'normal';
  end
end
