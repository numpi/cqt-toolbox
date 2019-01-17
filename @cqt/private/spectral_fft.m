function [um,up] = spectral_fft(vm,vp)
%SPECTRAL_FFT Computes the spectral factorization of a polynomial
%
%     [UM, UP] = SPECTRAL_FFT(VM, VP) computes the spectral factorization
%     of a polynomial P(Z) defined by the coefficients
%
%       P := [ VM(1) ... VM(end) VP(1) ... VP(END) ],
%
%     where the first element is the leading coefficient. That is, the
%     coefficient with the vectors reported above can be obtained as
%
%       P = conv(UM, UP)
%
%     The algorithm consists in evaluating the central 2m-1 coefficients of
%     the Laurent series 1/p(z) and then by computing the first and last
%     column of the inverse of the mxm symmetric Toeplitz matrix formed 
%     with these elements these columns, suitably scaled, provide an 
%     approximation to the desired factor.
%
%     The function does not use any special algorithm for solving the 
%     Toeplitz system.
%
%     June 2, 2016, By Dario A. Bini

% clean data
realflag = isreal(vm) && isreal(vp);

% vm = cln(vm); vp = cln(vp);
% compute the reciprocal of the Laurent polynomial
[tm,tp] = reciprocal_fft(vm,vp);
nm = length(vm); np = length(vp);  n = max(nm,np);
% Form the Toeplitz matrix
A = toep(tm,tp,n,n);
% Solve the Toeplitz system % it can be improved
b1 = zeros(n,1);  b2 = b1;  b1(1) = 1;  b2(n) = 1;
U = A\[b1,b2];
um = U(:,1);  up = U(end:-1:1,2);
% Clean vectors and normalize
% um = cln(um); up = cln(up);
n = min(length(um), length(up));
th = sum(um(1:n).*up(1:n));

if realflag
    ss = sign(vm(1) / th);
else
    ss = 1;
end

th = th * ss;

th = sqrt(vm(1)/th);
um = ss*th*um;  up = th*up;
end
