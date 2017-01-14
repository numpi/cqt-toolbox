function C = fqt_mldivide(A, B)
%FQT_MLDIVIDE

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt(0, 0, E.U, E.W(end:-1:1,end:-1:1), A.sz(1), size(E.U,2) + size(E.W,2));
E2 = cqt(0, 0, E.V, E.Z(end:-1:1,end:-1:1), A.sz(2), size(E.V,2) + size(E.Z,2));

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison
LUinv = inv(L) * inv(U);
LUB = LUinv * B;
LUE1 = LUinv * E1;

S = eye(size(E1, 2)) + full(E2.' * LUE1);

C = LUB - LUE1 * (S \ (E2.' * LUB));