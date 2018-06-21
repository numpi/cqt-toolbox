function C = qt_mrdivide(B, A)
%QT_MRDIVIDE Compute B * inv(A) for infinite CQT matrices

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
    C = inv(A) * B;
    return;
end

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, [], inf, size(E.U, 2));
E2 = cqt([], [], E.V, [], inf, size(E.V, 2));

% Compute A \ B by Sherman-Morrison

Linv = inv(L);
Uinv = inv(U);
LUinv = Linv * Uinv;
BLU = B * LUinv;
E2LU = E2.' * LUinv;

S = eye(size(E.U, 2)) + full(E2LU * E1);
E1SE2LU = -(E1 / S) * E2LU;
E1SE2LU.p = 1;
E1SE2LU.n = 1;
C = BLU * E1SE2LU;

end

