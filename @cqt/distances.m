function Dist = distances(A,x,N0,nsamp,advpx,dig)
% function Dist = distances(A,x,N0,nsamp,advpx,dig)
% This function computes the distances of each eigenvalue x(i) of the CQT matrix A
% from the closest eigenvalue of the truncated matrix A_N for N = 2^k N0, 
% k = 0,1,2,...,nsamp-1
% On Output
% The first line of the matrix Dist contains the truncation values N0,2*N0,...
% The remaining lines are such that Dist(i,j) is the  distance of
% eigenvalue i-1 of A from the closest value of A_N for N= N0*2^(j-1)

% D.A. Bini, March 2022

  [am, ap] = symbol(A);
  E = correction(A);

  if advpx
     mp.Digits(dig);
     am = mp(am); ap = mp(ap); E = mp(E); A=cqt(am,ap,E);
     tol = 10^(-dig+4);
  else
     tol =1.e-12;
  end
  
  nz = length(x);
  if advpx
     Dist = zeros(nz,nsamp,'mp');
  else
     Dist = zeros(nz,nsamp);
  end
  Dist(1,:) = N0*2.^(0:nsamp-1)';

  for j=1:nz
     N = N0;
     for k=1:nsamp
        A = sptoep(am,ap,sparse(E),N);
        if advpx
           B = mp(A-x(j)*speye(N));
        else
           B = A-x(j)*speye(N);
        end
        ei = invpower(B);
        Dist(j+1,k) = ei;
        N = N*2;
        if ei<tol
           break
        end
     end
  end 
end

function a = sptoep(am,ap,E,n)
% function a = toep(am,ap,E,n)
% The nxn matrix a = CQT(am,ap,E) is computed
  nm = length(am); np = length(ap);
  n = max(n,(nm+np)*2);
  a = am(1)*speye(n);
  for j=2:nm
    a = a+am(j)*diag(sparse(ones(n-j+1,1)),-j+1);
  end
  for j=2:np
    a = a+ap(j)*diag(sparse(ones(n-j+1,1)),j-1);
  end
  [k1,k2]=size(E);
  a(1:k1,1:k2) = a(1:k1,1:k2)+sparse(E);
end



function ei = invpower(A)
% function ei = invpower(A);
% This function provides the smallest modulus of the 
% eigenvalues of the sparse matrix A
% The inverse power method is used.
  maxit = 100;
  n = size(A,1);
  v=rand(n,1)-0.5+1i*(rand(n,1)-0.5);
  v = v/norm(v);
  sigold = 0;
  for j=1:maxit
     y = A\v;
     sig = sum(abs(y./v))/n;
     err = abs(sig-sigold)/sig;  
     if err<1.e-8  || sig>1.e-15
        ei = 1/sig;
        return
     end
     sigold = sig;
     v = y/norm(y);
  end
  fprintf('Warning: InvPower, Exceeded max number of iterations\n');
end


