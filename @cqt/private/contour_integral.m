function fM = contour_integral(f,x,r,max_it)
if ~exist('x', 'var')
    x=0;
end
if ~exist('r', 'var')
    r=1;
end
if ~exist('max_it', 'var')
    max_it = 64;
end
N = 4;
err = Inf;

threshold = sqrt(eps);
fM = f(x+r)*0; % Inizializza a zero tenendo conto dell'eventuale struttura cqt, h-matrice ecc
for j = 1 : N
    z = exp(2i * pi * j ./ (N));
    xi = r*z + x;
    fM = fM + (z * r / N) * f(xi);
end

while err > threshold && N <= max_it
    fnewM = .5 * fM;
    N = 2 * N;
    for j = 1 : 2 : N
        z = exp(2i * pi * j / N);
        xi = r*z + x;
        fnewM = fnewM + (z * r / N) * f(xi);
    end
    
    err = norm(fM - fnewM);
    fM = fnewM;
end