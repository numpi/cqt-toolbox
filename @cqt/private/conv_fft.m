function r = conv_fft(p, q)
%CONVFFT Compute conv(p, q) using the FFT.

%r = conv(p, q);
%return;

if ~isvector(p) || ~isvector(q)
	error('Arguments of CONV_FFT must be vectors');
end

p = reshape(p, 1, length(p));
q = reshape(q, 1, length(q));

l = length(p) + length(q) - 1;
p = [ p, zeros(1, l - length(p)) ];
q = [ q, zeros(1, l - length(q)) ];

r = ifft(fft(p) .* fft(q));

end

