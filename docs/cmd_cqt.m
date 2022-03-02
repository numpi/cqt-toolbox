%% Construction QT matrices
%
% The |cqt| command can be used to construct QT matrices from several input
% data. Different syntaxes are supported. 

%% Syntax
% * |T = cqt(neg, pos, E)| creates the semi-infinite cqt-matrix with
%      the symbol defined by the coefficients in |neg| and |pos|, which are 
%      expected to contain negative and positive powers of the Laurent polynomial, 
%      respectively. For consistency, it is checked that |neg(1) == pos(1)|. 
%      An optional finite correction |E| can be specified. 
% * |T = cqt(neg, pos, U, V)| creates the semi-infinite cqt-matrix with
%      the specified symbol (as above) and finite correction |U * V.'|. 
% * |T = cqt(neg, pos, A, B, m, n)| creates the $m \times  n$ QT matrix with
%   the specified symbol, top-left correction $A$ and bottom-right 
%   correction $B$. If one between |m| and |n| is |inf| then $B$ is ignored.
% * |T = cqt(neg, pos, U, V, W, Z, m, n)| creates the $m \times n$ 
%   QT matrix with the specified symbol and finite corrections |A = U * V.'|
%   and |B = W * Z.'|. If one between $m$ and $n$ is |inf| then $B$ is ignored.
% * |T = cqt(neg, pos)| creates the semi-infinite QT matrix with the 
%   specified symbol and an empty finite correction.
% * |T = cqt(A)| creates the semi-infinite cqt-matrix with a symbol 
%   equal to $0$ and finite correction equal to $A$. 
% * All the constructors can be called using the syntax:
%   |T = cqt('extended', arg1, ..., argk, v)|, 
%   which is equivalent to the commands
%   |T = cqt(arg1, ..., argk)| and |T = extend(T, v)|, 
%   which constructs a quasi-Toeplitz matrix plus a rank $1$ correction
%   of the form $ev^T$, where $e$ is the vector of all ones. 
% * |T = cqt('hankel', f)| constructs the Hankel matrix with symbol
%   given by |f|.

%% Example
%
A = cqt([ 1 2 3], [1 4 5 6], -ones(2));
A(1:8, 1:8)