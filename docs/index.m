%% User guide to the CQT toolbox
%
%% Introduction
%
% The CQT toolbox implements the arithmetic of $\mathcal{CQT}$ matrices, 
% which are defined as the set of matrices obtained summing a semi-infinite
% Toeplitz matrix and a correction to the top-left corner.
%
% A $\mathcal{CQT}$ matrix has the structure $A =  T(a(z))+ E_a$ where 
%
% $$ T(a(z)) := \left[ \begin{array}{cccc}
%            a_0 & a_1 & a_2 &  \cdots  \\
%            a_{-1} & a_0 &  \ddots & \ddots \\
%            a_{-2} &  \ddots & \ddots\\
%            \vdots & \ddots
%          \end{array} \right],
% $$
%
% and $E_a$ is a compact correction. Although $A$ is a matrix of infinite
% size, we assume that its symbol $a(z)$ admits a convergent Laurent series
% and that ${\rm vec}(E)$ is in $l^1(N)$. This allows to operate
% numerically on these objects by appropriately truncating these matrices.
%
% The CQT toolbox allows to represent and operate these matrices directly 
% in the MATLAB environment.
%
% If you are unfamiliar with QT matrices, visit the following link to
% understand how to interact with infinite matrices in MATLAB; later on,
% you may be interested in checking out some applications that can be
% treated in this framework. 
%
% * <QTDefinition.html First steps with QT matrices>
% * <Arithmetic.html Arithmetic operations between QT matrices>
% * <commands.html Overview of all commands and functions available in the toolbox>
% * <LinearSystems.html Solving linear systems with QT matrices>
%
%% Markov chains with infinite-dimensional state spaces
%
% QT matrices are a natural tool for analyzing Markov chains and random
% walks on infinite state spaces. We present a few examples of such
% analyses here.
% 
% * <JacksonExample.html Computing invariant probabilities> for certain 
%   Quasi-Birth-and-Death
%   Markov processes involving quadratic matrix equations.  
% * <QuadraticExample.html Solving a quadratix matrix equation> using 
%   cyclic reduction and CQT matrices. 
%
%%  Matrix functions of infinite matrices
%
% * <SquareRoot.html Computing the matrix square root> of semi-infinite
%   Toeplitz matrices with a low-rank correction. 
%
%% Computing eigenvalues of QT matrices
% * <doc_eig.html Computing eigenvalues> of QT matrices. 

