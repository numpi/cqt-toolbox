%% Distances from truncated matrices
%
% Given an infinite size QT matrix $A$, its finite truncation may (or may
% not) well approximate the spectrum of the infinite-dimensional operator. 
%
% The function |distance| computes the distance between the eigenvalues of
% the truncated matrix and the closest eigenvalue of the
% infinite-dimensional operator. 

%% Syntax
%
% |[D, x] = distances(A, algo, N0, nsamp, advpx, dig)|

%% Example
% See <doc_eig.html the documentation on computing eigenvalues> for further
% details on the use of this function. 