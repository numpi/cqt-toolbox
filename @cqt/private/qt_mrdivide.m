function C = qt_mrdivide(B, A)
%QT_MRDIVIDE Compute B * inv(A) for infinite CQT matrices

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices. 
E1 = cqt(0, 0, E.U);
E2 = cqt(0, 0, E.V);

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison
LUinv = inv(L) * inv(U);
BLU = B * LUinv;
E2LU = E2.' * LUinv;

S = eye(size(E.U, 2)) + correction(E2LU * E1);

C = BLU - BLU * (cqt(correction(E1) / S)) * E2LU;

end

