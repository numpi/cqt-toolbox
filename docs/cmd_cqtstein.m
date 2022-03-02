%% Solving Stein equations in QT format
%
% The theory for the solution of Stein equations is related with the one
% for solving Lyapunov and Sylvester equations, which can be done with the
% command <cmd_cqtlyap.m cqtlyap>. 
%
% A Stein equation has the form $AXB + X + C = 0$, for given matrices
% $A,B,C$. If the latter are $\mathcal{QT}$ matrices and the product of any
% two elements in the spectra of $A$ and $B$ is always smaller than $1 -
% \epsilon$ for some $\epsilon > 0$, then the solution $X$ is representable
% in $\mathcal{QT}$ format as well. 

%% Syntax
% * |X = cqtstein(A, B, C)| computes a matrix |X| that satisfies |A*X*B + X
%   + C = 0|, up to a small truncation error. 

%% Example
A = cqt([ .5 .2 ], [.5 .2 ]);
B = A;
C = cqt(1, 1);

X = cqtstein(A, B, C);

norm(A*X*B + X + C) / norm(C)