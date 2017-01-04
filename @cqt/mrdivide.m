function X = mrdivide(B, A)
%MRDIVIDE Compute B * inv(A). 

if isa(B, 'cqt') && ~isa(A, 'cqt')
    if isinf(B.sz(1))
        error('Incompatible dimensions');
    end
    
    X = cqt(0, 0, full(B) / A, [], B.sz(1), B.sz(2));
    
    return;
end

if isa(A,'cqt') && isa(B,'cqt')
    if isinf(A.sz(1))
        X = qt_mrdivide(B, A);
    else        
        X = fqt_mrdivide(B, A);
    end
else
	error('A and B must be cqt matrices');
end



end

