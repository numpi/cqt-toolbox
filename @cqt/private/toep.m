function C = toep(am,ap,n1,n2)
% function C = toep(am,ap,n1,n2)
% The function creates an n1xn1 Toeplitz matrix whose first column and row
% are am, and ap, respectively
a1 = length(am);  a2 = length(ap);
cm = zeros(n1,1); cp = zeros(n2,1);
cm(1:min(n1,a1))=am(1:min(n1,a1));
cp(1:min(n2,a2))=ap(1:min(n2,a2));
C = toeplitz(cm,cp);

