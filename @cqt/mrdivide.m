function X = mrdivide(B, A)
%MRDIVIDE Compute B * inv(A). 

if isa(B, 'cqt') && ~isa(A, 'cqt')
    if B.sz(2) ~= size(A,1)
        error('Incompatible dimensions');
    end
    rw = B.sz(1);
    if rw == inf
        B.sz(1) = max(size(B.U,1), length(B.n) + B.sz(2));
    end
    X = cqt([], [], full(B) / A, [], rw, B.sz(2));
    
    return;
end

if isa(A,'cqt') && isa(B,'cqt')
    if B.sz(2) ~= A.sz(1)
        error('Incompatible dimensions');
    end
    if isinf(A.sz(1))
        X = qt_mrdivide(B, A);
    else        
        X = fqt_mrdivide(B, A);
    end
else
	error('A and B must be cqt matrices');
end

end

