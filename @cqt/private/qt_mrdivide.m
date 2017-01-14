function C = qt_mrdivide(B, A)
%QT_MRDIVIDE Compute B * inv(A) for infinite CQT matrices

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, [], inf, size(E.U, 2));
E2 = cqt([], [], E.V, [], inf, size(E.V, 2));

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison
LUinv = inv(L) * inv(U);
BLU = B * LUinv;
E2LU = E2.' * LUinv;

S = eye(size(E.U, 2)) + full(E2LU * E1);

C = BLU - BLU * (E1 / S) * E2LU;

end

