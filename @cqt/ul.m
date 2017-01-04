function [U, L, E] = ul(A)
%UL Compute the UL factorization of A as A = U*L + E. 
%
% [U, L, E] = ul(A) computes the UL factorization of a CQT matrix A. The
% factors U, L, and E are all CQT matrices. U and L are upper and lower
% triangulat Toeplitz matrices, and E is a finite correction. 

switch cqtoption('inversion')
    case 'cr'
        spectral = @spectral_cr;
    case 'fft'
        spectral = @spectral_fft;
end

[l, u] = spectral(A.n, A.p);

L = cqt(l, l(1));
U = cqt(u(1), u);

% Set the correct sizes in L and U
L.sz = [ min(A.sz) , A.sz(2) ];
U.sz = [ A.sz(1) , min(A.sz) ];

% Compute the required correction E
if isinf(A.sz(1))
    [E1, E2] = correction(A);
    E = cqt(0, 0, E1, E2);
else
    [E1, E2, E3, E4] = correction(A);
    E = cqt(0, 0, E1, E2, E3, E4, A.sz(1), A.sz(2));
end

% In the finite case we need to add another term at the bottom of the
% matrix. 
if ~isinf(A.sz(1))
  % Compute the lower right corner correction
  hU = hankel(u(2:end));
  hL = hankel(l(2:end));
  
  m = min(length(u), length(l)) - 1;
  
  [ ~, new_u, new_l, ~ ] = lower_right_correction(u(1), u, l, l(1), ...
      A.sz(1) , min(A.sz), A.sz(2));  
  
  hU = hankel(new_u(2:end));
  hL = hankel(new_l(2:end));
  
  m = min(length(new_u), length(new_l)) - 1;
  
  hU = hU(:,1:m);
  hL = hL(1:m,:);  
  
  F = cqt(0, 0, 0, 0, hU(end:-1:1,end:-1:1), hL(end:-1:1,end:-1:1).', A.sz(1), A.sz(2));
  
  E = E + F;
end