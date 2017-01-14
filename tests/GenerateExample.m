function [T1, T2, U1, V1, sn1, sp1, U2, V2, sn2, sp2] = GenerateExample(n, k1, k2)
%GENERATEEXAMPLE Generate two random CQT matrices. 
%
% [T1, U1, V1, sn1, sp1, T2, U2, V2, sn2, sp2] = GENERATEEXAMPLE(n, k1, k2)
%     generates two CQT matrices T1 and T2 that can be used as test
%     examples. 
%
% Author: Leonardo Robol <leonardo.robol@cs.kuleuven.be>

[U1,~] = qr(randn(n, k1) + 1i * randn(n, k1)); [V1,~] = qr(randn(n, k1));
[U2,~] = qr(randn(n, k2) + 1i * randn(n, k2)); [V2,~] = qr(randn(n, k2));

sn1 = [ 15 rand(1, 5) ];
sp1 = [ sn1(1) , rand(1, 4) ];
sn2 = rand(1, 4);
sp2 = [ sn2(1), rand(1, 7) ];

T1 = cqt(sn1, sp1, U1, V1);
T2 = cqt(sn2, sp2, U2, V2);


end

