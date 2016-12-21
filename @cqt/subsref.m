function M = subsref(T, indices)
%SUBSREF Extract a submatrix of T, or access some properties. 

switch indices.type
    case '()'
        idx = indices.subs;
        if length(idx) ~= 2
            error('Specify two indices to slice a matrix');
        end
        
        row_indices = idx{1};
        col_indices = idx{2};
        
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
           end
           
        end
        
        
    otherwise
        error('Unsupported slicing operand');
end

end

