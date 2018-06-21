function M = slicemat(T, idx)
%SLICE Slice a portion of a CQT matrix.

if length(idx) ~= 2
    error('Specify two indices to slice a matrix');
end

row_indices = idx{1};
col_indices = idx{2};

if max(row_indices) > T.sz(1) || max(col_indices) > T.sz(2)
    error('Index exceeds matrix dimension')
end
if min(row_indices) <= 0 || min(col_indices) <= 0
    error('Negative indices are not supported')
end

if isconsecutive(row_indices) && isconsecutive(col_indices)
    M = slicemat_consecutive(T, ...
        row_indices(1), row_indices(end), ...
        col_indices(1), col_indices(end));
    return;
end

M = zeros(length(row_indices), length(col_indices));

for i = 1 : length(row_indices)
    for j = 1 : length(col_indices)
        if ( row_indices(i) > col_indices(j) ) && ...
                ( row_indices(i) - col_indices(j) < length(T.n) )
            M(i,j) = T.n(row_indices(i) - col_indices(j) + 1);
        elseif ( row_indices(i) <= col_indices(j) ) && ...
                ( col_indices(j) - row_indices(i) < length(T.p) )
            M(i,j) = T.p(col_indices(j) - row_indices(i) + 1);
        end
        
        if row_indices(i) <= size(T.U, 1) && ...
                col_indices(j) <= size(T.V, 1)
            M(i,j) = M(i,j) + T.U(row_indices(i), :) * T.V(col_indices(j), :).';
        end
        if row_indices(i) > T.sz(1) - size(T.W, 1) && ...
                col_indices(j) > T.sz(2) - size(T.Z, 1)
            M(i,j) = M(i,j) + T.W(T.sz(1) - row_indices(i) + 1, :) * T.Z(T.sz(2) - col_indices(j) + 1, :).';
        end
    end
    
end

if ~isempty(T.c)
    if length(T.c) < max(col_indices)
        T.c(max(col_indices)) = 0;
    end
    
    M = M + ones(size(M,1),1) * T.c(col_indices);
end

