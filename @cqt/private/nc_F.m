function [nc, G, exc] = nc_F(W, am, ap, G, x, variant,advpx)
  % function [nc, G, exc] = nc_F(W, am, ap, G, x, variant)
% Compute the Newton correction nc = det(A)/det'(A)
% for A = WV, if variant = 0     (Algorithm  2)
%     A = WV-xV if variant = 1   (Algorithm  4)
% where V = [I;G;G^2;...], G=F^p, and F is the Frobenius
% matrix associated with the polynomial a(z)-x
% Input
%    W,am,ap: parameters defining the problem
%    G: such that G=F^p, F Frobenius matrix associated with
%       the factor of a(z)-x with roots of modulus <1
%    x: initial approximation
%    variant: if 0, standard Newton correction, if 1, it applies a variant
%    advpx: if true, the toolbox Advanpix is used for high precision arithmetic
% Output
%    corr: correction such that y=x-corr is the new iteration
%    G: such that G=F^p, F Frobenius matrix associated with 
%       the factor of a(z)-y with roots of modulus <1
%    exc: true if exceeded the max number of cyclic reduction iterations

% By D.A. Bini, August 2, 2021

% Preliminaries

  exc = false;
  m = length(am)-1; 
  [h1,h2] = size(W);
  
  
% Compute p
  p = size(G,1);

% Enlarge W to a multiple of p
  cp = ceil(h2/p);
  N = cp*p; dif = N-h2;
  if dif>0
     if advpx
          W= [W,zeros(h1,dif,'mp')];
     else
          W= [W,zeros(h1,dif)];
     end
  end


% Compute the factor s and its derivative sp
  if advpx
    s = [-G(1,:).';mp('1')];
  else
    s = [-G(1,:).';1];
  end
  a = [am(m+1:-1:1);ap(2:end)]; 
  a(m+1)=a(m+1)-x;
  [sp,~,~] = deriv(m,a,s,advpx);
  if advpx
     sp=[sp;mp('0')];
  else
  sp=[sp;0];
  end

% compute the derivative of G
  Gp = derivG(s,sp,advpx);

% compute W(G) and its derivative, through G^k and (G^k)'
  A=W(:,1:p); Ap=W(:,p+1:2*p)*Gp;  Gip=Gp;
  if advpx
     Gi=eye(p,'mp');
  else
     Gi=eye(p);
  end



  for i=1:cp-2
     Gi = G*Gi; Gip = Gip*G+Gi*Gp;
     A = A+ W(:,i*p+1:i*p+p)*Gi;
     Ap = Ap+ W(:,(i+1)*p+1:(i+1)*p+p)*Gip;
  end
  A = A+W(:,(cp-1)*p+1:cp*p)*(G*Gi);
  
% compute the Newton correction



  den = trace(A\Ap);

  if variant ==1 
      den = den - m*sp(1)/s(1);
  end
  
  if isnan(den) 
     if advpx
       nc = mp('0');
     else
       nc = 0;
     end
  else
     if advpx
       nc = -mp('1')/den;
     else
       nc = -1/den;
     end
  end

% Compute G
  x = x-nc;
  p = m+wind(am,ap,x,advpx);
  % wi=wind(am,ap,x,advpx);
  [G, exc] = factorG(am,ap,p,x,advpx);

end

