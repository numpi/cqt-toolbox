function X = mldivide(A, B)

if isa(B, 'cqt') && ~isa(A, 'cqt')
    if B.sz(1) ~= size(A,2)
        error('Incompatible dimensions');
    end
    cl = B.sz(2);
    if cl == inf
	B.sz(2) = max(size(B.V,1), length(B.p) + B.sz(1));
    end
    X = cqt([], [], A \ full(B), [], size(A, 1), cl);
    
    return;
end

if isa(A,'cqt') && isa(B,'cqt')
    if isinf(A.sz(1))
        X = qt_mldivide(A, B); 
    else        
        X = fqt_mldivide(A, B);
    end
else
	error('A and B must be cqt matrices');
end
