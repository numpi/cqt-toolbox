function [U, S, V] = lanczos_svd(Afun, n, tol)
%LANCZOS_SVD Two-sided Lanczos for SVD.

V = randn(n, 1);
V(:,1) = V(:,1) / norm(V);

U = zeros(n, 0);

res = inf;

if ~exist('tol', 'var')
    tol = cqtoption('threshold');
end

alfa = [];
beta = [];

it = 0;

nrm = 0.0;
nrm_converged = false;

while res > tol
    it = it + 1;
    
    w = Afun(V(:,end), 'notrasp');
    
    % Reorthogonalize w against U
    w = w - U * (U' * w);
    w = w - U * (U' * w);
    
    alfa = [ alfa, norm(w) ];
	
	% keyboard
	
	if alfa(end) == 0
		res = 0;
		beta = [ beta, 0 ];
		break;
	end
    
    U = [ U , w / alfa(end) ];
    
    % Compute beta
    w = Afun(U(:,end), 'trasp');
    w = w - V * (V' * w);
    w = w - V * (V' * w);
    
    beta = [ beta, norm(w) ];
	
    V = [ V, w / beta(end) ];
	
	if size(V, 2) > size(V, 1)
		break;
	end	
    
    % Estimate the norm: if we have a good estimate, evaluate the
    % possibility of stopping the iteration.
    % fprintf(' Norm : %e, nrm_converged: %d, res = %e\n', beta(end) / nrm);
    if ~nrm_converged
        nrm_est = norm(diag(alfa) + diag(beta(1:end-1), 1));
        if (nrm_est - nrm) / nrm < cqtoption('threshold')
            nrm_converged = true;
        end
        nrm = nrm_est;
    else
        res = beta(end) / nrm;
    end
end

[Ul, S, Vl] = svd(diag(alfa) + diag(beta(1:end-1), 1));

% Possibly perform recompressione
rk = sum(diag(S) > tol * S(1,1));

if rk == 0
	U = [];
	S = [];
	V = [];
else
	U = U * Ul(:,1:rk);
	V = V(:,1:end-1) * Vl(:,1:rk);
	S = S(1:rk, 1:rk);
end

end

