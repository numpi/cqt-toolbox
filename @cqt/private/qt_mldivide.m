function C = qt_mldivide(A, B)
%QT_MLDIVIDE Compute inv(A) * B. 

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices. 
E1 = cqt(0, 0, E.U);
E2 = cqt(0, 0, E.V);

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison
LUinv = inv(L) * inv(U);
LUB = LUinv * B;
LUE1 = LUinv * E1;

S = eye(size(E.U, 2)) + correction(E2.' * LUE1);

C = LUB - LUE1 * cqt(S \ correction(E2.' * LUB));
