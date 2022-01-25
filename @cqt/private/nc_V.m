function [nc,r] = nc_V(W, am, ap, r, x, variant,advpx)
% function [nc,r] = nc_V(W, am, ap, r, x, variant)
% Compute the Newton correction nc det(A)/det'(A)
% Input:
%    W,am,ap parameters defining the problem
%    r roots of a(z)-x
%    x initial approximation
%    variant: either 0 for A = WV
%             or 1 for A = HV-x V
%             where V is the Vandermonde matrix associated with r 
%    advpx: if true, the toolbox Advanpix is used for high precision arithmetic

% Output
%    corr: correction such that y=x-corr is the new iteration
%    r: roots of a(z)-y

% By D.A. Bini, August 2, 2021

% Parameters
  
  
  np = length(ap);
  m = length(am)-1;  n = np-1;
  h2 = size(W,2);
  p = length(r);

% compute the derivative xp of x
if advpx
  xp = r*mp('0');
else
    xp = r*0;
end
  for i =1:n
     xp = xp + i*ap(i+1)*r.^(i-1);
  end
  for i =1:m
     xp = xp - i*am(i+1)*r.^(-i-1);
  end
  if advpx
    xp = mp('1')./xp;
  else
    xp = 1./xp;
  end


% compute the matrix V
  if advpx
     V = ones(h2,p,'mp');
  else
     V = ones(h2,p);
  end
  for j=2:h2
     V(j,:) = V(j-1,:).*r;
  end

% compute the derivative of V  
  Vp = V; Vp(2:end,:)=V(1:end-1,:);  
  if advpx
     Vp(1,:)=mp('0'); 
  else
     Vp(1,:)=0; 
  end
if advpx
  D = diag(mp(0:h2-1));
else
  D = diag(0:h2-1);
end
  Vp = D*Vp;
  D = diag(xp); Vp = Vp*D;
  A = W*V; Ap = W*Vp;

% Newton correction
  den = trace(A\Ap);
  if variant==1
     den = den - m*sum(xp./r);
  end
  if isnan(den)
     nc = 0;
  else
     if advpx
        nc = mp('1')/den;
     else
        nc = 1/den;
     end
  end

% Update r
  x = x-nc;
% compute the roots
  a = [ap(end:-1:1);am(2:end)];
  b = a; b(np) = b(np)-x;
  ro = roots(b);
  [s1,s2] = sort(abs(ro));
  r = ro(s2).';
% select roots in the unit disk
  ind = find(s1<1);
  p = length(ind);
  r = r(1:p);
end

