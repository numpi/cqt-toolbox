%% Solving Jackson waiting Queues. 
%
%% Generating the problem
%
% We use the function CQTGALLERY to construct the problems from the paper
% of Motyer and Taylor [1]. In this case we choose some of the examples. 

i = 3;

[Am1, A0, A1, hA0] = cqtgallery('jackson', i);

%% 
% We will need a certain threshold to regulate the truncation of the
% output, which can be set here. 

threshold = 1e-12;

%% 
% We solve the associated matrix equations $A_{-1} + A_0 G + A_1 G^2 = 0$ 
% and $R^2  A_{-1} + R A_0 + A_1 = 0$ by means of cyclic reduction. 

[G, R] = cr(Am1, A0, A1, 10);

% Let's double check if these are the correct solutions
fprintf('Residue of the right solution G: %e\n', norm(Am1 + A0 * G + A1* G^2));
fprintf('Residue of the left solution R: %e\n', norm(R^2 * Am1 + R * A0 + A1));

%% 
% We are now ready to compute the first $m$ entries of the invariant
% probability vector $\pi$, which we store in $\pi_0$. We choose $m$ a
% priori, but the parameter needs to be adjusted so that the accuracy is
% high enough. 

M = hA0 + A1 * G;
[Mm, Mp] = symbol(M);
m = 4 * max([ size(correction(M)), length(Mm), length(Mp) ]);

% Get a vector in the kernel
[Q, ~, ~] = qr(M(1:m,1:m));
pi0 = Q(:,end); 

% Check if we computed the vector accurately enough (otherwise you might
% want to increase m). 
res = cqt(pi0') * M;
fprintf ('Residual norm on pi0: %e\n', norm(res));

%% 
% We can now compute the other rows of the vector $\pi$ by multiplying
% $\pi_0$ on the right by the matrix $R$. Notice that these are all
% infinite vectors, therefore we will need to store them truncated in the
% output. 

infinite_pi0 = cqt(pi0');
pi = pi0;

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

fprintf('Residue on the computed pi: %e\n', norm(R));
