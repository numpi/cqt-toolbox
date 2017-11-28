function [U, S, V] = random_svd(Afun, n, tol)
%RANDOM_SVD Compute an approximate SVD with random sampling. 

if ~exist('tol', 'var')
    tol = eps;
end

% Oversampling parameter
p = 10;

% Estimated rank 
k = 5;

converged = false;

while ~ converged
    [Q, ~] = qr(Afun(randn(n, k + p), 'notrasp'), 0);
    B = Afun(Q, 'trasp')';
    
    [U, S, V] = svd(B, 'econ');
    
    if S(end,end) / S(1,1) < tol
        converged = true;
    else
        k = 2 * k;
    end
end

rk = sum(diag(S) > tol * S(1,1));

U = Q * U(:,1:rk);
S = S(1:rk,1:rk);
V = V(:,1:rk);

end

