function C = qt_mldivide(A, B)
%QT_MLDIVIDE Compute inv(A) * B.

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
	C = inv(A) * B;
	return;
end

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, [], inf, size(E.U, 2));
E2 = cqt([], [], E.V, [], inf, size(E.V, 2));

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison

%LUinv = inv(L) * inv(U);
%LUB = LUinv * B;
%LUE1 = LUinv * E1;
Linv = inv(L);
Uinv = inv(U);
LUB  = Linv * (Uinv * B);
LUE1 = Linv * (Uinv * E1);

S = eye(size(E.U, 2)) + full(E2.' * LUE1);

C = LUB - LUE1 * (S \ E2.') * LUB;
