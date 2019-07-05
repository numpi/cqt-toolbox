function [U, V] = hankel_compress(a, b, strategy)
%HANKEL_COMPRESS Compress the product of two Hankel matrices.

n = max(length(a), length(b));

a = [ a , zeros(1, n - length(a)) ];

if ~isempty(b)
	b = [ b , zeros(1, n - length(b)) ];
end

if ~exist('strategy', 'var')
    strategy = 'random';
end

switch strategy
    
    case 'lanczos'
		if isempty(b)
			[U, S, V] = lanczos_svd(@(x,t) hankel_matvec1(a, x), n);
		else
			[U, S, V] = lanczos_svd(@(x, t) hankel_matvec(a,b,x,t), n);
		end
        
        U = U * sqrt(S);
        V = V * sqrt(S);
        
    case 'random'
		if isempty(b)
			[U, S, V] = random_svd(@(x,t) hankel_matvec1(a, x), n);
		else
			[U, S, V] = random_svd(@(x, t) hankel_matvec(a,b,x,t), n);
		end
        
        U = U * sqrt(S);
        V = V * sqrt(S);
        
    otherwise
        error('Unsupported compression strategy selected');
end

end

function y = hankel_matvec1(a, x)
	h = length(a);
	y = toepmult_fft(a(end), a(end:-1:1), h, h, x(end:-1:1));
end

function y = hankel_matvec(a, b, x, trasp)
%HANKEL_MATVEC Perform the matrix multiplication y = A * B * x

if strcmp(trasp, 'trasp')
    y = hankel_matvec(b, a, x, 'notrasp');
else
    h = min(length(a), length(b));
    
    % Perform the multiplication y = B * x
    y = toepmult_fft(b(h:-1:1), b(h), h, h, x);
    
    % ... and then y = A * (B * x);
    y = toepmult_fft(a(h), a(h:-1:1), h, h, y);
end
end

