function M = subsref(T, indices)
%SUBSREF Extract a submatrix of T, or access some properties.

switch indices.type
    case '()'
        M = slicemat(T, indices.subs);
    otherwise
        error('Unsupported slicing operand');
end

end

