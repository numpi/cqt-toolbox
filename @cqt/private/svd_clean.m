function [U, V] = svd_clean(U, V, S, nrm)

epsilon = nrm * cqtoption('threshold');

if ~isempty(U) || ~isempty(V)
    [U, V, ~] = svd_clean_rec(U, V, S, epsilon);
end

end

function [U, V, S, epsilon] = svd_clean_rec(U,V, S, epsilon)

    if isempty(U) || isempty(V)
        U = zeros(0, length(S));
        V = zeros(0, length(S));
        return;
    end
    
    nrmU = norm(U(end,:) .* S);
    nrmV = norm(V(end,:) .* S);

    if min(nrmU, nrmV) < epsilon
        if nrmU < nrmV
            epsilon = epsilon - nrmU;
            U = U(1:end-1,:);
        else
            epsilon = epsilon - nrmV;
            V = V(1:end-1,:);
        end
        
        [U, V, S, epsilon] = svd_clean_rec(U, V, S, epsilon);
    end
end
