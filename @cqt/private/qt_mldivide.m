function C = qt_mldivide(A, B)
%QT_MLDIVIDE Compute inv(A) * B.

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
    C = inv(A) \ B;
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
LUB  = LUinv * B;
LUE1 = LUinv * E1;

S = eye(size(E.U, 2)) + full(E2.' * LUE1);

LUE1SE2 = - LUE1 * (S \ E2.');
LUE1SE2.p = 1;
LUE1SE2.n = 1;
C = LUE1SE2 * LUB;
