function [f, u, ff] = fact(b,n,verbose)
  % b has n zeros inside the unit disk and m outside
  % b(z) = b_0+b_1z+...+b_Nz^N
  % f(z) = f_0+f_1z+...+f_{n-1}z^{n-1}+z^n
  % u(z) = u_0+u_1z+...+u_mz^m
  % N = n+m
  % f is stored in an n-dimensional vector
  % u is stored in an (m+1)-dimensional vector

% Bruno Iannazzo Jan 2022
  ff = {};
  N = length(b)-1;
  if size(b,1)==1
    b = b.';
  end
  b0 = b(1:n);
  b1 = b(n+1:N+1);
  m = N-n;

  if n>m
    error('This works only for n<=m')
  end

  % initial value
  col = flipud(b0);
  row = [0 b1'];
  
  F2 = toeplitz([col;zeros(10,1)],[row zeros(1,10)]);
  F2 = F2(1:(n+3),1:(n+3));  %det(F2),
  F3 = toeplitz(b,zeros(1,n));
  f1 = F2\[1;zeros(n+2,1)]; f1 = f1(1:n);
  f = F3(1:n,1:n)*f1;
  
  
  %F2 = toeplitz(col,row);
  %F2 = F2(1:n,1:n);
  %F3 = toeplitz(b,zeros(1,n));
  %f = F3(1:n,1:n)*(F2\[1;zeros(n-1,1)]);
  %ff{1}=f;

  % iteration
  maxiter = 200;
  for ell = 1 : maxiter
    % phi
    F0 = [tril(toeplitz(f)) zeros(n,m+1-n)];
    F1 = [triu(toeplitz([1 fliplr(f') zeros(1,m-n)]))];
    u = F1\b1;
    phi = b0 - F0*u;
    % phi'
    u0 = u(1:n);
    u1 = u(n+1:m+1);
    G0 = F0(:,1:m-n+1);
    G1 = F1(1:m-n+1,1:m-n+1);
    d = u0 - G0*(G1\u1);
    phi1(:,1) = d;
    Z = eye(n); Z=Z(:,2:n);
    F = [Z, -f];
    for h = 2 : n
      d = F*d;
      phi1(:,h) = d;
    end
    % step
    cor = phi1\phi; 
    err = norm(cor./f,inf);
    f = f + cor;
   % ff{ell+1} = f;
    if verbose
        fprintf('W-H fact, Newton: ite=%d err=%e\n',ell, err);
    end
    if err<1.e-15
       return
    end
    % criterion
  end

