%% User guide to the CQT toolbox
%
%% The basics
%
% The CQT toolbox implements the arithmetic of $\mathcal{CQT}$ matrices, which are
% defined as the set of matrices obtained summing a semi-infinite Toeplitz
% matrix and a correction to the top-left corner.
%
% More precisely, we say that $X$ is $\mathcal{CQT}$ if $X = T + E$ with
% $T$ semi-infinite Toeplitz and $E$ is a matrix such that ${\rm vec}(E)$
% is in $l^1(N)$.
%
% The CQT toolbox allows to represent these matrices directly in MATLAB.
%
% * <ArithmeticExample.html First steps with CQT matrices>
%
%% Some more involved examples
% 
% * <JacksonExample.html Computing invariant probabilities> for certain Quasi-Birth-and-Death
%   Markov processes involving quadratic matrix equations.  
% * <SquareRoot.html Computing the matrix square root> of semi-infinite
%   Toeplitz matrices with a low-rank correction. 
% * <QuadraticExample.html Solving a quadratix matrix equation> using cyclic
%   reduction and CQT matrices. 