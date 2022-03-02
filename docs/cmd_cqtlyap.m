%% Lyapunov and Sylvester equations
%
% If $A,B,C$ are $\mathcal{QT}$ matrices such that the spectra of $A$ and
% $-B$ are separated, then the solution of the matrix equation $AX + XB +
% C= 0$ is again a $\mathcal{QT}$ matrix. 
%
% The solution can be computed numerically by calling the function
% |cqtlyap|. 

%% Syntax
% |X = cqtlyap(A, B, C)| returns |X| that satisfies |A*X+X*B+C=0|, up to a
%   small truncation error. 

%% Example
%
% We generate a matrix $A$ which is positive definite, so that setting $B =
% A$ guarantees the separation of the spectra between $A$ and $-B$. We
% solve with $C = I$, the identity. 

A = cqt([ 2.5, 1 ], [ 2.5, 1 ], rand(4));
B = A;
C = cqt(1, 1); 

X = cqtlyap(A, B, C);

norm(A*X + X*B + C) / norm(C)