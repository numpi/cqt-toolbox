function [U, V] = svd_clean(U, V, S, nrm)

epsilon = nrm * cqtoption('threshold');

if ~isempty(U) || ~isempty(V)
    [U, V, ~] = svd_clean_loop(U, V, S, epsilon);
end

end

% This function contains two versions of svd_clean: a recursive one, as
% documented in the paper "Quasi-Toeplitz matrix arithmetic: A MATLAB
% toolbox", and a loop-based one, which is equivalent but more efficient,
% and is used by default. 
%
% The old function is left for historical and documentation purposes. 

function [U, V, S] = svd_clean_loop(U, V, S, epsilon)

    acc_error = 0.0;
    
    % Scale U and V to make everything simpler later on
    UU = U * diag(S);
    VV = V * diag(S);
    
    % Default step size: this will be reduced when we get close to finding
    % the right truncation level
    step_size = 32;
    
    if isempty(U) || isempty(V)
        U = zeros(0, length(S));
        V = zeros(0, length(S));
        return;
    end
    
    u_ptr = size(U, 1);
    v_ptr = size(V, 1);
    
    while acc_error < epsilon
        u_start = max(1, u_ptr - step_size + 1);
        v_start = max(1, v_ptr - step_size + 1);
        
        % Estimating the 2-norm by 1-norm // not sharp, but much faster
        nrmU = norm(UU(u_start:u_ptr,:), 1) * sqrt(step_size);
        nrmV = norm(VV(v_start:v_ptr,:), 1) * sqrt(step_size);
        
        if min([ nrmU, nrmV ]) > (epsilon - acc_error)
            if step_size == 1
                % No further truncation is possible, exist from here
                U = U(1:u_ptr, :); V = V(1:v_ptr, :); return;
            else
                step_size = step_size / 2;
                continue;
            end
        end
        
        if nrmU < nrmV
            u_ptr = u_ptr - step_size;
            acc_error = acc_error + nrmU;
        else
            v_ptr = v_ptr - step_size;
            acc_error = acc_error + nrmV;
        end
        
        if min([ u_ptr, v_ptr ]) == 0
            k = length(S);
            U = zeros(0, k); V = zeros(0, k); return;
        end
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
