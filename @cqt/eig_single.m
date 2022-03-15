function [x1, U, y, info] = eig_single(A, x0, varargin)
%function [x1, U, y, info] = eig_single(A, x0, varargin)
% Compute an eigenvalue of A = cqt(am,ap,E) by means of
% Newton iteration starting from an initial guess x0
% Input:
%   A = cqt(am,ap,E)
%   x0 : initial guess
% Optional values:
% 'algo', a,  where a = 1,2,3,4, according to the algorithm used:
%       1 Newton's iteration applied to det(WV)=0 (Vandermonde version)
%       2 Newton's iteration applied to det(W[I;G;G^2;...)=0 (Frobenius version)
%       3 Newton's iteration applied to det(HV-x V)=0 (Vandermonde version)
%       4 Newton's iteration applied to det(H(G)-x I)=0 (Frobenius version)
%       default value: 1
% 'maxit', m, where m is the maximum number of iterations, default m=20
% 'epsilon', ep, where ep is the relative precision for the halt criterion,
%        default ep=1.e3*eps
% 'verbose', ver, if ver=true, some information is printed at run time,
%        default ver=false
% 'residual', res, if res=true, the residual error is computed, default res=false
% 'advpx', ad, if ad=true, the multiprecision toolbox advanpix is used, 
%        default ad=false
% 'digits', dig, where dig is the number of decimal digits precision used in 
%        advanpix, default  dig=34
% 'eigenvector', k, where k is the number of block components of the eigenvector
%        default value k=0
%
% Output:
% x1 : eigenvalue
% if 'eigenvector', 0, then
% U : is either the pxp Vandermonde or G=F^p, where F is the Frobenius matrix
% y : is the vector such that [VDy; VD^2y; VD^3y; ...] or [Gy; G^2y; G^3y;...], 
%     is the corresponding eigenvector
% if 'eigenvector', k, where  k > 0, then
% U : is the eigenvector and y is empty
% info : this variable contains further information, more specifically,
% info.flg : 0 No values in the connected component containing x1 are eigenvalues
%              (the number p of zeros of a(z)-x1 with modulus <1 is p=0)
%            1 p=q x1 is an isolated eigenvalue  (q is the number of equations)
%            2 p>q All values in the connected component containing x1  are eigvals
%            3 p<q x1 is an isolated eigval
%            4 p<q x1 is not eigenvalue
%            5 p<q x1 is outside the connected component containing x0
%            6 diverging sequence
%            7 exceeded the max number of iterations
% info.res :  residual error in the first q components of the equation Av-x1 v=0
% info.vres : vector with the residual errors per step 
% info.vresest : vector with the residual errors per step estimated by the SVD 
% info.it : number of iterations


% By D.A. Bini, L. Robol, January 28, 2022

info = struct;
p = inputParser;
addOptional(p, 'algo', 1);
addOptional(p, 'maxit', 20);
addOptional(p, 'epsilon', []);
addOptional(p, 'verbose', false);
addOptional(p, 'residual', false);
addOptional(p, 'advpx', false);
addOptional(p, 'digits', 34);

% Number of components of the eigenvector to compute. 0 means to avoid
% computing it and return U, y instead. 
addOptional(p, 'eigenvector', 0);

parse(p, varargin{:});
algo = p.Results.algo;
maxit = p.Results.maxit;
epsil = p.Results.epsilon;
verbose = p.Results.verbose;
residual = p.Results.residual;
advpx = p.Results.advpx;
digits = p.Results.digits;
eigenvector = p.Results.eigenvector;

if isempty(epsil)
    if advpx
        epsil = mp(10)^(-digits+3);
    else
        epsil = 1e3 * eps;
    end
end

% Parameters
  warning off; 		% warning on;
  [am, ap] = symbol(A);
  E = correction(A);  
  if advpx
      am = mp(am); ap = mp(ap); E = mp(E);
      x0 = mp(x0);
  end

% adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end
  if am(1)~=ap(1)
      fprintf('WARNING: am(1) and ap(1) have inconsistent values\n')
      return
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
  am = am(1:nm); ap = ap(1:np); A = cqt(am,ap,E);
  
  if m==0
    fprintf('The matrix is block upper triangular, use eig_all instead\n');
    x1 = [];
    return
  end
  if n==0
    fprintf('The matrix is block lower triangular, use eig_all instead\n');
    x1 = [];
    return
  end

 if advpx
    mp.Digits(digits);
    am = mp(am);
    ap = mp(ap);
    E = mp(E);
    A = cqt(am,ap,E);
 end
 
% create the polynomial a(z) and check sizes   % March 7 2022
  a = [ap(end:-1:1);am(2:end)];
  np = length(ap);
  nm = length(am);  m = nm-1;
  [h1,h2] = size(E);  
  if h1<m
     if advpx
         E=[E;zeros(m-h1,h2,'mp')];
     else
         E=[E;zeros(m-h1,h2)];
     end
     h1 = m;
  end
% compute an upper  bound to the moduli of the eigenvalues
  norm1 = norm(a,1);
  mxval = norm1 + norm(E,inf);

% compute the rank of the exceeding part of E
  if h1>m
     E2 = E(m+1:end,:);
     r2 = rank(E2);
     [~,R,P] = qr(E2);
     R = R*P';
     R = R(1:r2,:);
  else
     if advpx
        r2 = mp('0');   R=mp('0');
     else
        r2 = 0;   R=0;
     end
  end
  q = m + r2;

% create the triangular Toeplitz matrix B and the matrix W
  if advpx
     e1 = zeros(m,1,'mp'); 
     W = zeros(q,m+h2,'mp'); 
  else
     e1 = zeros(m,1); 
     W = zeros(q,m+h2);
  end 
  e1(1) = am(m+1); 
  if length(e1)>1    %%%%%%%%%%%%%%% March 7 2021
     B = toeplitz(e1,am(m+1:-1:2));
  else
     B = e1;
  end
  W(1:m,1:m) = -B;
  W(1:m,m+1:m+h2) = E(1:m,:); 
  if r2>0
      W(m+1:m+r2,m+1:end) = R;
  end  

% Preliminaries for the iterations
  x1 = []; info.flg = []; info.res = []; info.vres = []; info.vresest=[];
  it = []; U = []; Y = []; y = []; r = [];
    
  if mod(algo,2)==1	% Vandermonde version
%    compute the roots
     b = a; b(np) = b(np)-x0;
     ro = roots(b);
     [s1,s2] = sort(abs(ro));
     r = ro(s2).';
%    select roots in the unit disk
     ind = find(s1<1);
     p = length(ind);
     r = r(1:p);
  else 	        % Frobenius version
%    Compute p and enlarge W and Rw to a multiple of p
     w = wind(am,ap,x0,advpx); 
     p = w + m;
     hm2 = h2+m;
     cp = ceil(hm2/p);
     N =cp*p; dif = N-hm2;
     if dif >0
        if advpx
            W = [W,zeros(q,dif,'mp')];
        else
            W = [W,zeros(q,dif)];
        end
     end
%    Compute G
     G = factorG(am,ap,p,x0,advpx);
  end     

% Start the iterations 
  errold=1.e300;
  info.vres = zeros(maxit,1); info.vresest = info.vres;
  for it = 1:maxit 
     if verbose
        fprintf('step = %d, q = %d,  p= %d, x0=%e %+ei\n', it, q, p, real(x0),imag(x0));
     end 

%-1  distinguish the cases p=0, p>q, p=q, p<q
     if p==0		   % x is not eigenvalue
        info.flg = 0;  x1 = x0;
        if verbose
           fprintf('return 1, flg=0: No eigenvalues in this component\n');
           fprintf('No eigenvalues in the component containing x0 =%d+%di\n', real(x0), imag(x0));
        end
        info.it = it;
        return % nov 21
     elseif p>q	   % Continuous set of eigenvalues
        info.flg = 2;  x1 = x0;
	if verbose
           fprintf('return 2,  flg=2: x belongs to a continuous set of eigenvalues\n');
           fprintf('x0 is eigval together with all the values in this connected component\n');
           x1 = x0;
	end
        info.it = it;
        return
     elseif p<q		% x might be (unilkely) eigenvalue
        if verbose
           fprintf('case     p<q\n');
        end
        S = W(1:p,:);     % Nov 21
        info.flg = 3;
     else		% p=q
        S = W; % Nov 21 
        info.flg = 1;
     end

%-2  select the algorithm 
     switch algo	
        case 1
           [corr,r] = nc_V(S,am,ap,r,x0,0,advpx);
        case 2
           [corr,G] = nc_F(S,am,ap,G,x0,0,advpx); 
        case 3
           [corr,r] = nc_V(S,am,ap,r,x0,1,advpx); 
        case 4
           [corr,G] = nc_F(S,am,ap,G,x0,1,advpx); 
        otherwise
           fprintf('Wrong algorithm selection\n')
           if verbose
		fprintf('return 3: Wrong algorithm selection\n');
	   end
           info.it = it;
           return
     end    
     x1 = x0-corr;   
     err = abs(corr)/abs(x0);
     x0 = x1;
     if verbose
        fprintf('      rel_correction = %1.5e, flag = %d\n', err, info.flg);
     end

%-3  Update p
     if mod(algo,2)==1    	% Vandermone form
         pnew = length(r);
     else                   % Frobenius form
         pnew =size(G,1);
     end
     if pnew ~= p
         if verbose
             fprintf('return 4 pnew>p: change of component\n');
         end
         info.flg = 5;
         info.res =1.e300;
         info.it = it;
         return
     end

%-4  Update U and prepare to compute y and the residual error
     if mod(algo,2)==1    	% Vandermone form
%        compute the Vandermonde matrix U
         mh2 = m+h2;  
         if advpx
            U = ones(mh2,p,'mp');
         else
            U = ones(mh2,p);
         end
         for j=2:mh2
            U(j,:) = U(j-1,:).*r;
         end
     else			% Frobenius form
         mh2 = m+h2;
         cp = ceil(mh2/p);
         if advpx
            U = zeros(p*cp,p,'mp');
            U(1:p,1:p) = eye(p,'mp');
         else
            U = zeros(p*cp,p);
            U(1:p,1:p) = eye(p);
         end  
         for j=1:cp-1 
             U(p*j+1:p*j+p,:) = U(p*(j-1)+1:p*j,:)*G;
         end
     end
     
        
%    Compute the vector y by means of SVD, and the residual error
   
     Y = W*U;                % Nov 21  
     [~,sy,vy]=svd(Y); sy = diag(sy);
     y=vy(:,end);            % (y is beta)     
     resalt = sy(end)/sy(1); % nov 27 2021   alternate residual  
     if residual
        info.res = qtresidual(A,x1,U,y,algo,advpx);
        info.vres(it) = info.res;         % store the residual at each step
     end
     info.vresest(it) = resalt;   % store the alternate residual
     if verbose
        if residual
           fprintf('      residual = %1.5e, resalt = %1.5e \n', info.res, resalt);
        else
           fprintf('      resalt = %1.5e \n',  resalt);
        end
     end

% Check convergence: stp=true if convergence occurr

     stp = err < epsil;
     stp1 = (err < 1000*epsil) && errold<=err;
     stp = stp || stp1;
     errold = err;
     
     if stp   % Check flag before halting
        if info.flg == 1 || info.flg == 3
            % Check if the eigenvector is requested
            if eigenvector > 0
              U = qteigenvector(eigenvector, U, y, length(am) - 1, algo, advpx);
              y = [];
            end
        end

        if info.flg == 1    % x is eigenvalue with p=q
            if verbose
               fprintf('return 4, flg=1: x is eigenvalue with p=q\n');
            end
            info.it = it;
            return            
        elseif info.flg == 3      % overdetermined problem
            if info.res<epsil
               if verbose
                  fprintf('return 6, flg=3: x is eigenvalue and p<q\n');
               end
               info.it = it;
               return
            else
               info.flg = 4;
               if verbose
                  fprintf('return 7, flg=4: x is not eigenvalue and p<q\n');
               end
               info.it = it;
               return
            end
        end
     end 
     if abs(x0)>mxval*10^2  % Check if diverging sequence
        x1 = 0;info.flg = 6; 
	    if verbose
           fprintf('return 7, flg=6: diverging sequence\n');
	    end
        info.it = it;
        return
     end
  end  % end of iterations
  
  if verbose 
       fprintf('warning: flg=7, reached the max number of iterations, res = %e\n',info.res)
       fprintf('x=%e %+ei\n',real(x0),imag(x0))
  end
  x1 = x0; info.flg = 7; info.it = it;
  



