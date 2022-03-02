%% Limit of the columns
%
% For Extended QT matrices, it is of interest to compute the limit for $i
% \to \infty$ of $A_{ij}$, for a fixed $j$. 
%
% While for standard QT matrices $A = T(a) + E$ this is always zero, for
% matrices that also have a correction of the form $ev^T$ this is equal to
% the vector v. 

%% Syntax
% * |v = limit(A)|

%% Example
%
% To obtain some interesting data, we need to use extended QT matrices. 

A = cqt('extended', [1 2], [1 2], 1, [1 2 3 4 5 6]);
limit(A)