function G = factorG(am,ap,p,x,advpx)
% function G = factorG(am,ap,p,x,advpx)
% Compute G = F^p where F is the Frobenius matrix associated with the
% roots of the Laurent polynomial a(z) of modulus less than 1
% by means of solving a matrix equation using the U-based method
% advpx: if true, the toolbox Advanpix is used for high precision computation

% By D.A. Bini, December, 2021 

%  nwtwh = strcmp(cqtoption('wiener-hopf'), 'newton');
   nwtwh = false;
 % adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end
  % apply cr-qbd 
  a = [am(end:-1:1);ap(2:end)];
  nm = length(am);
  a(nm) = a(nm)-x;
  if nwtwh
     [g,~] = newton_WH(a,p,advpx);
 %    U = triu(toeplitz(g(1:p))); L = tril(toeplitz(g(p+1:-1:2)));
 %    G = -L\U;  % Barnett identity
     F = diag(ones(p-1,1),1);
     F(p,:) = -g(1:p).';
     G = F^p;
    
     
  else
     G = mycr(a,p,advpx);
  end
    
 

