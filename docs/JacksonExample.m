%% Solving Jackson waiting Queues.
%
%% Generating the problem
%
% We use the function CQTGALLERY to construct the problems from the paper
% of Motyer and Taylor [1]. In this case we choose the example with index
% equal to 2.  We also start a timer using the function TIC, so we will
% measure how long it takes to solve the entire problem. 

tic;

[Am1, A0, A1, hA0] = cqtgallery('jackson', 2);

%%
% We will need a certain threshold to regulate the truncation of the
% output, which can be set here.

threshold = 1e-12;

%%
% We solve the associated matrix equations $A_{-1} + A_0 G + A_1 G^2 = 0$
% and $R^2  A_{-1} + R A_0 + A_1 = 0$ by means of cyclic reduction.

[G, R] = cr(Am1, A0, A1, 12);

% Let's double check if these are the correct solutions
fprintf('Residue of the right solution G: %1.2e\n', norm(Am1 + A0 * G + A1* G^2));
fprintf('Residue of the left solution R: %1.2e\n', norm(R^2 * Am1 + R * A0 + A1));

%%
% We can now compute the matrix $M = A_0 + A_1 G$ so that it is block
% tridiagonal and block Toeplitz, with at most the top-left block different
% from the other diagonal ones. We need to find an integer $m$ such that
% the correction is contained in the first $m$ rows and columns, and such
% that the numerical band of $M$ is less than $m$.

M = hA0 + A1 * G;
[Mm, Mp] = symbol(M);

m = max([ size(correction(M)), length(Mm), length(Mp) ]);

%%
% Construct the blocks $\hat M_0, M_{-1}, M_0, M_1$ so that
%
% $$ M = \left[ \begin{array}{cccc} \hat M_0 & M_1 && \\ M_{-1} & M_0 & M_1
% & \\ & \ddots & \ddots & \ddots  \end{array} \right] $$

hM0 = M(1:m,1:m);
M1  = M(1:m,m+1:2*m);
M0  = M(m+1:2*m,m+1:2*m);
Mm1 = M(m+1:2*m,1:m);

%%
% Compute the small $G_M$ and $R_M$ relative to the block Toeplitz system
% with blocks $M_i$.

[GM, RM] = cr(Mm1, M0, M1);

%%
% Compute a vector in the left kernel of the matrix $\hat M_0 + M_1 G$, and
% construct all the other components of the vector $\pi_0$ which are not
% negligible.

[Q, ~] = qr((hM0 + M1 * GM));
pi0 = Q(:,end)';

pim1 = pi0;
nrm_pi0 = norm(pi0, inf);

while norm(pim1, inf) > nrm_pi0 * eps
    pim1 = pim1 * RM;
    pi0 = [ pi0, pim1 ];
end

% Check if we computed the vector accurately enough (otherwise you might
% want to increase m).
res = cqt(pi0) * M;
fprintf ('Residual norm on pi0: %e\n', norm(res));

%%
% We can now compute the other rows of the vector $\pi$ by multiplying
% $\pi_0$ on the right by the matrix $R$. Notice that these are all
% infinite vectors, therefore we will need to store them truncated in the
% output.

% Normalization
s = cqt(pi0) / (cqt(1,1) - R);
s = sum(s(1,1:1024));
pi0 = pi0 / s;

infinite_pi0 = cqt(pi0);
pi = pi0';

while norm(infinite_pi0) > threshold
    infinite_pi0 = infinite_pi0 * R;
    pi = [ pi , infinite_pi0(1, 1:size(pi,1))' ];
end

%%
% We can now check the residual on the computed vector.

E = cqt([], [], 1);
J = cqt(0, [0 1]);

pi = cqt(pi)';
R = E * pi * (hA0 - A0) + J * pi * Am1 ...
    + J' * pi * A1 + pi * A0;

fprintf('Residue on the computed pi: %1.2e\n', norm(R));
fprintf('Residue on the computed pi (infty norm): %1.2e\n', norm(R(1:4096, 1:4096), inf));

total_time = toc;
fprintf('Total time used: %4.2f s\n', total_time);