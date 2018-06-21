%% Computing the matrix square root
%
% We consider the computation of the matrix square root described 
% in the paper
%
% "Quasi-Toeplitz matrix arithmetic: a MATLAB toolbox", by D. A. Bini, S.
% Massei, L. Robol, Numerical Algorithms, 2018.
%
% We construct a semi-infinite Toeplitz matrix $A$ with symbol
%
% $$ a(z) = \frac{1}{4} \cdot \left( z^{-2} + z^{-1} + 4 + 2z + 
% z^2 \right),  $$
%
% and we want to compute $B = A^{\frac{1}{2}}$. Let us start by defining
% the symbols 
am = [1 .25 .25];
ap = [1 .5 .25];

%%
% We then choose a correction with a support of 1024 entries, but only rank
% 3. 

corr_size = 1024; 

U = randn(corr_size, 3);
V = randn(corr_size, 3);

%%
% Normalize to avoid making the eigenvalues too close to the origin

U = U / norm(U * V') / 5;
    
A = cqt(am, ap, U, V);
tic; B = sqrtm(A); toc; 

%%
% And finally, we check if the result is correct (up to some threshold). 
norm(A - B^2) / norm(A)