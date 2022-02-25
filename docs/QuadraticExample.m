%% Solving a quadratic matrix equation
%
% QT matrices are particularly useful in the numerical treatment of random
% walks over the positive orthant ${\bf N}^2$. Suppose that a walker
% can move over all points of non-negative integer coordinates, and:
%
% * moves to the left or right with probabilities $l$ and $r$;
% * moves to the top and bottom with probabilities $t$ and $b$; 
% * stays in the current entry with probability $s$;
% * is given suitable boundary conditions when reaches $\{ 0 \} \times \bf N$ or 
%   ${\bf N} \times \{ 0 \}$.
%
% If we consider a semi-infinite matrix $X_t$ that represents the probability 
% of being in a given state of indices $(i,j)$ at a given time $t$, then after
% one discrete time step the walker will have be in a new state with
% probabilities
%
% $$ X_{t+1} = sX_t + 
%  \left[\begin{array}{cccc} 
%   & b && \\ t && b & \\ & t & & \ddots \\ & & \ddots &  \\
%   \end{array}\right] X_t
%    + X_t \left[\begin{array}{cccc} 
%   & r && \\ l && r & \\ & l & & \ddots \\ & & \ddots &  \\ \end{array}\right]
%   = sI + VX + XH, 
%  $$
%  where $V$ encodes vertical movements, and $H$ the horizontal ones. 
%    
% Under these assumptions, the probability transition matrix for the Markov
% chain can be written as follows, up to rearranging the entries in $X_t$ by stacking
% a sequence of (infinite) vectors one of top of the other:
%
% $$
%   P = sI + H \otimes I + I \otimes V = \left[ \begin{array}{cccc}
%     \hat A_0 & A_1 && \\ 
%      A_{-1} & A_0 & A_1 & \\ 
%      & A_{-1} & \ddots & \ddots \\
%      & & \ddots & \ddots \\
%   \end{array} \right], \qquad 
%   A_0 = sI + V, \ A_1 = r\cdot I, \ A_2 = l\cdot I.
% $$
%
% In particular, the blocks $A_i$ are infinite tridiagonal QT matrices, and $P$ has the
% same block structure. The invariant probability vector of the Markov
% chain can be computed by the Matrix Analytic Method by M. Neuts [1],
% which requires to solve the matrix equation $ X = A_{-1} + A_0 + A_1 X^2
% $.
%
% We consider the matrix equation associated with a QBD process and a
% random walk on $\{1 \ldots \infty\} \times \{ 1 \ldots \infty\}$. 
%
% We now show how to obtain an equation of this form, that may be obtained
% by allowing even more freedom enabling diagonal movements (which gives a
% tridiagonal structure in $A_1$ and $A_{-1}$ as well).
%
% To begin with, we define three CQT matrices with symbols 
%
% $$ a_{1}(z) = \frac{1}{4} \left( 2z^{-1} + 2 + z \right), \qquad 
%    a_{0}(z) = \frac{1}{10} \left(2z^{-1} + z \right), \qquad 
%    a_{-1}(z) = \frac{1}{6} \left( 6z^{-1} + 3 + 4z \right). $$

a1n = [ .5 .5 ];
a1p = [ .5 .25];

a0p = [ 0 0.1 ];
a0n = [ 0 0.2 ];

am1n = [ .5 1 ];
am1p = [ .5 2/3 ];

%% 
% Then, we need to scale these matrices in order to enforce
% row-stochasticity of the matrix polynomial $\varphi(z) = z^{-1} A_1 + A_0
% + z A_1$ when evaluated at $z = 1$. 

rowsum = sum([ a1n, a1p(2), a0n, a0p(2),  am1n, am1p(2) ]);
a1n = a1n / rowsum; a1p = a1p / rowsum;
a0n = a0n / rowsum; a0p = a0p / rowsum;
am1n = am1n / rowsum; am1p = am1p / rowsum;

%%
% Now that we know the symbols, we can create the CQT objects. 
%
A0 = cqt(a0n, a0p, a0n(2), a0p(2));
A1 = cqt(a1n, a1p, a1n(2), a1p(2));
Am1 = cqt(am1n, am1p, am1n(2), am1p(2));

%%
% We now solve the quadratic equation $A_{-1} + A_0 G + A_1 G^2 = G$ using
% cyclic reduction (which can be called using the CR function in the
% toolbox), and we verify the solution. 
%
% Note that the function CR actually solves $A_{-1} + A_0 G + A_1 G^2 = 0$,
% so we need to shift the middle coefficient with the identity. 

tic; G = cr(Am1, A0 - cqt(1, 1), A1); toc
norm(Am1 + A0 *G + A1 * G^2 - G)

%% Using Newton iteration
% The same problem can be solved using the Newton iteration, which is
% implemented in the function QUADNEWT. 
tic; G = quadnewt(Am1, A0 - cqt(1, 1), A1); toc
norm(Am1 + A0 *G + A1 * G^2 - G)
