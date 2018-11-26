function [ym,yp] = reciprocal_fft(am,ap)
% function [ym,ya] = reciprocal_fft(a)
% compute the coefficients of the reciprocal of the Laurent polynomial
% am(1)+am(2)/z+am(3)/z^2+... +ap(2)z +ap(2)z^2+... by using evaluation/
% interpolation at the roots of 1. The number of interpolation points is 
% set adaptively.
% On output, ym, yp are the coefficients of the negative/positive  powers,
% ym(1)=yp(1) is the constant coefficient
% June 22, 2016. By Dario A. Bini

a1 = length(am);  a2 = length(ap);

am = reshape(am, 1, a1);
ap = reshape(ap, 1, a2);

% constants:
maxiter = 18;
epsi = 1.e-15;
realflag = isreal(am) && isreal(ap);

% compute the first approximation
k = 1 + ceil(log(a1 + a2)/log(2));
Na = 2^k; N = Na;
% fill the vector to size N
v = zeros(1,N);
cntr = N/2+1;
v(cntr:-1:cntr-a1+1) = am;  v(cntr+1:cntr+a2-1) = ap(2:end);
% compute the scaling diag for permutation
scl = 2*mod([1:N],2)-1;
% evaluation-interpolation
s = 1./fft(scl.*v);
y = ifft(s);   y = y./scl;
if realflag
    y = real(y);
end
ym = y(cntr:-1:1);   yp = y(cntr:end);
% compute the next approximations until convergence
for iter=1:maxiter
    k = k+1;
    N = 2*N;
    % fill the vector to size N
    v = zeros(1,N);
    cntr = N/2+1;
    v(cntr:-1:cntr-a1+1) = am;  v(cntr+1:cntr+a2-1) = ap(2:end);
    % update the scaling diag for permutation
    scl = [scl,scl];
    % evaluation-interpolation
    s = 1./fft(scl.*v);
    yn = ifft(s); yn = yn./scl;
    
    if realflag
        yn = real(yn);
    end
    ynm = yn(cntr:-1:1);   ynp = yn(cntr:end);
    % Check convergence
    erm = norm(ym - ynm(1:length(ym)),'inf');
    erp = norm(yp - ynp(1:length(yp)),'inf');
    ym = ynm;  yp = ynp;
    if erm <epsi*norm(ym,'inf') && erp<epsi*norm(yp,'inf')
        % Clean the output
        % ym = cln(ym); yp = cln(yp);
        break
    end
end

if realflag
    ym = real(ym);
    yp = real(yp);
end

if iter>=maxiter
    disp([ 'Warning: reciprocal_fft has reached the max number' ...
        'of iterations. The error is' ])
    disp([erm,erp])
    %pause
end

