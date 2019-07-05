function y = hankel_matvec(a, x)
%HANKEL_MATVEC Given a Hankel matrix H with symbol A, compute Y = H * X. 

n = size(x, 1);

y = toepmult_fft(a(end), a(end:-1:1), n, n, x(end:-1:1, :));

end

