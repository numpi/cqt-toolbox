function [s,u]=newton_WH(a,p,advpx)
% function [s,u]=newton_WH(a,p,advpx)
% compute the Wiener-Hopf factorization of a=[a0,...,an]
% by means of Newton iteration starting from s equal to
% the polynomial having the same roots of a in the unit disk
% and random u
% p is the number of roots in the unit disk
epsi = eps*10^2; 

if size(a,1)==1
   a = a.';
end

if advpx
  epsi = 10^(3-mp.Digits);
end
n=length(a)-1; ph=n-p;
w=zeros(n+1,1);
eu1=zeros(p,1);eu1(1)=1;
es1=zeros(ph+1,1);es1(1)=1;
if advpx
  w=zeros(n+1,1,'mp');
  eu1=zeros(p,1,'mp'); eu1(1)=mp(1);
  es1=zeros(ph+1,1,'mp'); es1(1)=mp(1);
end
% Choosing the starting approximation
r = roots(a(end:-1:1));
ind = find(abs(r)<1);
nr =length(ind);
if p~=nr
   fprintf('Warning in Newton_HR: mismatching in the value of p, p=%d, nr=%e\n',p,nr);
   p = nr;
end
s = [-r(ind(1)),1].';
for i=2:nr
  s = conv(s,[-r(ind(i)),1].');
end

u = rand(ph+1,1)-0.5; u(end)=a(end);
if advpx
   s = mp(s); u = mp(u);
end
x = [s(1:end-1);u];
for it = 1:20
% Sylvester matrix
  w1 = w;  w1(1:ph+1) = u;
  M = toeplitz(w1,eu1*u(1));
  w1 = w;
  w1(1:p+1) = s; 
  N = toeplitz(w1,es1*s(1));
  Sylv = [M,N];
  c = conv(s,u)-a;
  corr = Sylv\c;  % Newton correction
  x = x-corr;
%  err=norm(corr,inf)/norm(x,inf);
  err = norm(corr./x,inf);
  s = x (1:p); s(p+1) = 1;
  if advpx
      s(p+1) = mp(1);
  end
  u = x(p+1:end);
  if err<epsi
     return
  end
end
