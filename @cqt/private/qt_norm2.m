function nrm = qt_norm2(A)
%QT_NORM2 Norm 2 estimator

[am, ap] = symbol(A);
n = max([ length(am), length(ap), size(correction(A)) ]);

v = randn(n, 1);

At = A';

s = 0;

for i = 1 : 30
    olds = s;
    s = norm(v);
    
    if abs(olds - s) < abs(s) * 1e-5
        break;
    end
    
    v = v / s;
    w = v;
    w = A * cqt(w);
    w = correction(At * w);
    v = w;
end

nrm = sqrt(s);

end

