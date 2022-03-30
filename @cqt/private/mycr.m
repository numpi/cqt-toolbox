function [G, exc] = mycr(a,p,advpx)
% function [G, exc] = mycr(a,p,advpx)
% Applies Cyclic Reduction to compute the minimal solution G
% of the equation Am1 + A0 X + A1 x^2 +...+ AM x^{M-1} = 0
% Am1=toeplitz([a_0,0,...,0],[a_0,....a_{p-1}],
% A_0=toeplitz([a_p,...,a_1],[a_p,...a_{2p-1}]);
% A_1=toeplitz([a_{2p},...,a_{p+1}],[a_{2p+1},...a_{3p-1}])
% the most convenient size partitioning kxk is chosen, over the set
% k = j+max(p,length(a)-p-1), j=0,1,....,2*p
% advpx: if true, the Advanpix toolbox is used for high precision arithmetic
% In output: exc is true if the maximum number of iterations has been reached

% Dario Bini March 6, 2021.
  exc = false;
  N = length(a); debug = true; debug = false;
  epsi = 1.e-16; max_it = 30;
  if advpx
     epsi = 10^(-mp.Digits);
  end
  k = max(p+1,N-p+1);
% find the size with the best cond
  am = double(a(p+1:-1:1));
  ap = double(a(p+1:end));
  bm =zeros(k+p,1);
  bp=bm;
  bm(1:p+1)=am;
  bp(1:length(ap))=ap; 
  T = toeplitz(bm,bp);
  c = zeros(p,1);
  for i=0:p
     ni = k+i-1;
     A= T(1:ni,1:ni);  %toeplitz(bm(1:ni),bp(1:ni));
     c(i+1)=1/rcond(A); %cond(A);
     if c(i+1)<1.e2
       break
     end
  end
  c=c(1:i+1);
  mc = min(c); im = find(c==mc);
  k = k+im(1)-1;
  
  %%%%%%%%%%%%%%%% 6/3/22
  k=k-1;
  %%%%%%%%%%%%%%%
  
  if advpx
     aa = zeros(3*k,1,'mp');
     e1 = zeros(k,1,'mp');
  else
     aa = zeros(3*k,1);
     e1 = zeros(k,1);
  end     
  aa(k-p+1:k-p+N)=a;
% Form the matrices Am1, A0, A1
  e1(1) = aa(1);
  S = toeplitz(e1,aa);
  Am1 = S(:,1:k); 
  A0 = S(:,k+1:k+k); 
  A1 = S(:,2*k+1:3*k);
 
% Start CR
  hA0 = A0; AAm = Am1;
  
  for it=1:max_it
     AA0 = A0 \ [ A1, Am1 ];
     temp1 = AA0(:, 1:k);
     temp2 = AA0(:, k+1:end);
     temp3 = Am1 * temp1;
     temp4 = A1 * temp2;   
     hA0 = hA0 - temp4;
     A0 = A0 - temp4 - temp3;
     Am1 = -Am1 * temp2;
     A1 = -A1 * temp1;   
     nrm_m1 = norm(Am1, inf);
     nrm_1  = norm(A1, inf) ;
     err = sqrt(nrm_m1*nrm_1);
     if debug
        fprintf('CR: step=%d, err=%e\n',it,err);
     end
     if err < epsi 
        break
     end   
  end
  if it==max_it
      disp('CR: maximum number of iterations reached!')
      exc = true;
  end
  G = - hA0 \ AAm;
  G = G(1:p,end-p+1:end);
 

