%% Solving a quadratic matrix equation
%
% We consider the matrix equation associated with a QBD process and a
% random walk on $\{1 \ldots n\} \times \{ 1 \ldots \infty\}$. 
% These are finite matrices ($n  \times n$). 
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
% Now that we know the symbols, we can create the CQT objects. We create
% them finite of dimension $n$, and in this case we choose $n = 2^{15}$. 
n = 2^15;
A0 = cqt(a0n, a0p, a0n(2), a0p(2), n, n);
A1 = cqt(a1n, a1p, a1n(2), a1p(2), n, n);
Am1 = cqt(am1n, am1p, am1n(2), am1p(2), n, n);

%%
% We now solve the quadratic equation $A_{-1} + A_0 G + A_1 G^2 = G$ using
% cyclic reduction (which can be called using the CR function in the
% toolbox), and we verify the solution. 
%
% Note that the function CR actually solves $A_{-1} + A_0 G + A_1 G^2 = 0$,
% so we need to shift the middle coefficient with the identity. 

tic; G = cr(Am1, A0 - cqt(1, 1, [], [], n, n), A1); toc
norm(Am1 + A0 *G + A1 * G^2 - G)

%% Using Newton iteration
% The same problem can be solved using the Newton iteration, which is
% implemented in the function QUADNEWT. 
tic; G = quadnewt(Am1, A0 - cqt(1, 1, [], [], n, n), A1); toc
norm(Am1 + A0 *G + A1 * G^2 - G)
