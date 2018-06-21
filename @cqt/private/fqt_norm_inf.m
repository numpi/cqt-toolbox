function r = fqt_norm_inf(T)
%FQT_NORM_1 Compute the 1-norm of T.

r = 0;

m = size(T, 1);
n = size(T, 2);

% Compute the row with index size(U,1) + 1, if necessary
if size(T.U,1) + size(T.W,1) < m
    i = size(T.U,1) + 1;
    i1 = min(i, length(T.n));
    i2 = min(n-i+1, length(T.p));
    
    ll = sum(abs([ T.n(1:i1) , T.p(2:i2) ]));
    r = max(r, ll);
    
    for i = size(T.U,1) + 2: m - size(T.W,1)
        if i > length(T.n)
            break;
        end
        
        ll = ll + abs(T.n(i));
        if n -i + 2 < length(T.p)
            ll = ll - abs(T.p(n-i+2));
        end
        r = max(r, ll);
    end
    
    if size(T.U,1) + size(T.W,1) < n
        r = max(r, norm([ T.n, T.p(2:end) ], 1));
    end
end

for i = 1 : min(size(T.U,1), m - size(T.W,1))
    indices = [ 1 : size(T.V,1), ...
        max(size(T.V,1)+1, i - length(T.n) + 1) : min(n, i + length(T.p)-1) ];
    r = max(r, norm(slicemat(T, { i, indices }), 1));
end

for i = m - size(T.W, 1) + 1 : size(T.U, 1)
    indices = [ 1 : size(T.V, 1), ...
        max(size(T.V,1)+1, i - length(T.n) + 1) : min(n - size(T.Z,1), i + length(T.p)-1), ...
        (n - size(T.Z,1) + 1) : n ];
    r = max(r, norm(slicemat(T, { i,indices }), 1));
end

for i = m - size(T.W,1) + 1 : m
    indices = [ max(1, i - length(T.n) + 1) : min(n - size(T.Z,1), i + length(T.p)-1), ...
        (n - size(T.Z,1) + 1) : n ];
    r = max(r, norm(slicemat(T, { i, indices }), 1));
end


end

