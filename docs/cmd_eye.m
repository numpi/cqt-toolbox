%% Constructing and identity in QT format
%
% The command |eye| constructs an identity represented in the QT format. 

%% Syntax
% * |I = eye(inf, 'like', cqt)|
% * |I = eye(n, 'like', cqt)|

%% Example
%
% We can use the command to build the identity for both finite and infinite
% matrices, as follows:

I = eye(inf, 'like', cqt);
I(1:4, 1:4)