function r = fqt_norm_2(T)
%FQT_NORM_2 Estimate the 2-norm of T by the power method. 

m = size(T, 1);
n = size(T, 2);

x = randn(n, 1);
x = x / norm(x);

diff = inf;
r = 0;

threshold = 1e-6;

while log(abs(diff)) > threshold
    y = full(T' * (T * cqt(0, 0, x, 0, n, 1)));
    
    oldr = r;
    r = norm(y);
    
    if oldr == 0
        diff = inf;
    else
        diff = r / oldr;
    end
    
    x = y / r;
end

r = sqrt(r);

end

