%% Diagonal extraction from QT matrices
%
% Given a QT matrix $A$ with finite or infinite size, the command |diag| 
% extracts another QT matrix that only contains diagonal entries. Thanks to
% the Toeplitz structure, this vector will be almost constant, with the 
% only exception of the top and bottom entries due to the compact correction. 

%% Syntax
% * |d = diag(A)|

%% Example
%
% We build a $6 \times 6$ QT matrix with corrections only in the entries in
% position $(1,1)$ and $(6,6)$, extract the diagonal part and convert it to
% a full matrix using the <cmd_full.html full> command. 

A = cqt([1 2], [1 2], -1, -2, 6, 6);
d = full(diag(A))