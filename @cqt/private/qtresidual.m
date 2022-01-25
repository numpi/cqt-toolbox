function res = qtresidual(A,x,V,y,algo,advpx)
% function res = qtresidual(am,ap,E,x,V,y,algo,advpx)
% Compute the residual error ||Av-lambda v||
% Input
%   The QT matrix A
%   The parameters x,V,y,r,algo provided by the function qteigen
%   advpx: if true, the toolbox Advanpix is used for high precision arithmetic
% Output
%   res: the residual error ||Av-lambda v||/(||v||*||A||)

% By D.A. Bini, August 2021

  [am, ap] = symbol(A);
  E = correction(A);

  p = length(y);
  h1 = size(E,1); h2 = size(E,2);
  m = length(am)-1; n = length(ap)-1;
  n1 = max(m,h1); n2 = max(n1+n,h2);
  N = ceil(n2/p);
  v = qteigenvector(N,V,y,m,algo,advpx);
  v = v(1:n2);
  if advpx
     H = zeros(n1,n2,'mp');
     am1=zeros(n1,1,'mp'); ap1=zeros(n2,1,'mp');
  else
     H = zeros(n1,n2);
     am1=zeros(n1,1); ap1=zeros(n2,1);
  end
  H(1:h1,1:h2) = E;
  if h1<=m
      k = m;
  else
      k = m+1;
  end
  am1(1:k) = am(1:k);  ap1(1:n+1) = ap;
  am1(1) = am1(1)-x;  ap1(1) = ap1(1)-x;
  T = toeplitz(am1,ap1);
  H(:,1:size(T,2)) = H(:,1:size(T,2))+T;
  w = H*v;
  res = norm(w)/norm(v);
  res = res/(norm(am,inf)+norm(ap(2:end),inf)+norm(E,inf)); %%% Ago 2021
  return

