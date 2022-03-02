%% Evaluation of generic matrix functions
%
% The function |funm| evaluates a generic matrix function $f(A)$ for some
% $f(z)$ analytic over the ball centered at zero and of radius equal to the
% norm of $A$. 
%
% The approach uses a contour integral on the boundary of such ball to
% evaluate the function. 

%% Syntax
% * |X = funm(A, f)| evaluate the function $f(z)$ at $A$. 
% * |X = funm(A, f, 'opt1', value1, 'opt2', value2, ...| allows to specify
% further options. In particular, |'max_it'| is the maximum number of
% quadrature points to use, and |'poles'| is a list of the poles of the
% function in case it's meromorphic; these are taken into account by
% computing the required correction to the function. 

%% Example
%
% As an example we evaluate the matrix exponential, and compare with the
% (much more efficient) implementation in |expm|. 

A = cqt([ 2 1 1], [ 2 1 ]);
eA = funm(A, @(z) exp(z));
eA2 = expm(A);

res = norm(eA - eA2) / norm(eA2)