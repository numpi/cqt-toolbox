%% Matrix square root
%
% The function |sqrtm| evaluates the matrix square root of a $\mathcal{QT}$
% matrix $A$. This can be done using either the Denman-Beaver iteration, or
% the Cyclic Reduction. In most cases, the two algorithms should deliver
% similar results, although the Cyclic Reduction may be preferable for
% matrices closer to singularity. 

%% Syntax
% * |X = sqrtm(A)| evaluate the matrix square root. 
% * |cqtoption('sqrt', method')| selects the method used for the
% computation (see <cmd_cqtoption.html cqtoption>); the valid options are
% either |'cr'| or |'db'|. 

%% Example
%
% In the first example we take a fairly well-conditioned matrix, and rely
% on Denman-Beavers (the default algorithm).
A = cqt([ 2 .8 ], [ 2 .8 ]);
sA = sqrtm(A);
res = norm(sA^2 - A) / norm(A)

%%
% We may consider a matrix closer to singularity, and compare the accuracy
% of Denman-Beavers vs the Cyclic Reduction.
A = cqt([ 2 .99 ], [ 2 .99 ]);
cqtoption('sqrt', 'db');
sA = sqrtm(A);
resDB = norm(sA^2 - A) / norm(A)

cqtoption('sqrt', 'cr');
sA = sqrtm(A);
resCR = norm(sA^2 - A) / norm(A)

