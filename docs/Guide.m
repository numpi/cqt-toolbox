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
% <ArithmeticExample.html First steps with CQT matrices>
%
% <JacksonExample.html Solving QBD processes>