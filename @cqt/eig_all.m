% function [x, xcont, res, it, ei] = eig_all(AA, algo, maxit, epsi, fact, verbose, plotfig, residual, advpx, digits)
function [x, xcont, res, it, ei] = eig_all(AA, varargin)

% function [x1, xcont, res, it, ei] = eig_all(A, algo, maxit, epsi, fact, verbose, plotfig, residual, advpx, digits)
% Compute the eigenvalue of the QT matrix A by functional iteration
% starting from the eigenvalues of a finite section of A of suitable size 
% In Input:
% am, ap, E parameters defining the QT matrix A
% algo: 1--8 algorithm number
% fact: factor for determining the size of the finite section n =fact*max(input size)
% epsi: error bound
% In output:
% x: vector with the isolated eigenvalues  (p=q)
% xsusp: suspect eigenvalues, i.e. maxit exceeded and res<1.e-8
% xcont: eigenvalue in a continuous set
% res: vector with the residual errors of x
% cnd: vector with the condition number of x
% it: vector with the number of iterations
% ei: vector with the eigenvalues of the finite section

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
        epsi = 10^(-digits+3);
    else
        epsi = 1e3 * eps;
    end
end

% By D.A. Bini, December 27, 2021

[am, ap] = symbol(AA);
E = correction(AA);

% adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end

  n = max([length(am), length(ap),size(E,1),size(E,2)]);
  n = n*fact;
  A = toep(single(am),single(ap),single(E),n);
  tic;  ei = eig(A);  teig = toc;
  ei = double(ei);
  if advpx
     mp.Digits(digits);
     ei = mp(ei);
  end

  tic;
  c1 = 1; c2 = 1; 
  x = []; xcont = []; res=[]; it=[];
  if am(1)~=ap(1)
      fprintf('WARNING: am(1) and ap(1) have inconsistent values\n')
      return
  end

  for j=1:n
     w = wind(am,ap,ei(j),advpx);
     if w>=0  %%%%%%%%      
         [z,~,~,info]=eig(AA, ei(j), 'algo', algo, 'maxit', maxit, ...
                'epsilon', epsi*100, 'verbose', verbose, ...
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
  if  ~isempty(x)      % refine the approximations   
    maxit = 2; epsix=1.e-20;
    if verbose
       fprintf('\n qteig_all: refining the approximations....\n');
    end
    for j=1:length(x)
          [x(j),~,~,info]=eig(AA, x(j), 'algo', algo, 'maxit', maxit, ...
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
     qtrange(am,ap);
     drawnow;
     hold off;
  end
  fprintf('Time required by eig = %e time required by qteig =  %e\n',teig,tqt); 
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


function a = toep(am,ap,E,n)
% function a = toep(am,ap,E,n)
% The nxn matrix a = CQT(am,ap,E) is computed
  nm = length(am); np = length(ap);
  n = max(n,(nm+np)*2);
  a = am(1)*eye(n);
  for j=2:nm
     a = a+am(j)*diag(ones(n-j+1,1),-j+1);
  end
  for j=2:np
     a = a+ap(j)*diag(ones(n-j+1,1),j-1);
  end
 [k1,k2] = size(E);
  a(1:k1,1:k2) = a(1:k1,1:k2)+E;
end

