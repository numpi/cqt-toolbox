function X = mldivide(A, B)

if isa(B, 'cqt') && ~isa(A, 'cqt')
    if isinf(B.sz(1))
        error('Incompatible dimensions');
    end
    
    X = cqt(0, 0, A \ full(B), [], size(A, 1), B.sz(2));
    
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
