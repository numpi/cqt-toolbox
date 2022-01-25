function [sd,ud,u] = deriv(m,a,s,advpx)
% function [sd,ud,u] = deriv(m,a,s,advpx)
% Given a=[a0,a1,...,aN]-lambda[0,...,1,0,...0] and the factor [s0,s1,...,sm]
% compute u such that a = conv(u,s)
% compute the first derivative sd and ud of s and u respectively
% with respect to lambda
% advpx: if true the toolbox Advanpix is used for high precision computation

% By D.A. Bini, March 3, 2021

% compute u 
  na = length(a); ns = length(s);
  nu = na -ns+1;
  n = 2^ceil(log2(na));
  if advpx
    at = [a;zeros(n-na,1,'mp')];
    st = [s;zeros(n-ns,1,'mp')];
  else
    at = [a;zeros(n-na,1)];
    st = [s;zeros(n-ns,1)];
  end

  fa = fft(at);
  fs = fft(st);
  u = ifft(fa./fs);
  u = u(1:nu);

% compute sd, ud
% Form matrices U and S, then the Sylvester resultant R
  M = nu+ns-1;
  
  if advpx
    eu1 = zeros(ns,1,'mp');eu1(1)=u(1);
    es1 = zeros(nu,1,'mp');es1(1)=s(1);
    u1 = [u;zeros(M-nu,1,'mp')];
    s1 = [s;zeros(M-ns,1,'mp')];
  else
    eu1 = zeros(ns,1);eu1(1)=u(1);
    es1 = zeros(nu,1);es1(1)=s(1);
    u1 = [u;zeros(M-nu,1)];
    s1 = [s;zeros(M-ns,1)];
  end
  U = toeplitz(u1,eu1);
  S = toeplitz(s1,es1);
  
% restrict the matrices S and U
  U = U(1:end-1,1:end-1);
  S = S(1:end-1,1:end-1);  
  
% form the Sylvester resultant matrix R
  R = [U,S];
  
% solve the system
  b = zeros(M-1,1);
  if advpx
    b(m+1)=-mp('1'); % nov 21
%    b(ns)=-mp('1');
  else
    b(m+1)=-1;       % nov 21
%    b(ns)=-1';
  end
  
  der = R\b;
  sd = der(1:ns-1);
  ud = der(ns:end);
