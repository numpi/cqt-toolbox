function C = fqt_mrdivide(B, A)
%FQT_MRDIVIDE Compute B * inv(A) for finite CQT matrices

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt(0, 0, E.U, E.W(end:-1:1,end:-1:1), A.sz(1), size(E.U,2) + size(E.W,2));
E2 = cqt(0, 0, E.V, E.Z(end:-1:1,end:-1:1), A.sz(2), size(E.V,2) + size(E.Z,2));

LUinv = inv(L) * inv(U);

% Compute A \ B by Sherman-Morrison
LUinv = inv(L) * inv(U);
BLU = B * LUinv;
E2LU = E2.' * LUinv;

S = eye(size(E1, 2)) + full(E2LU * E1);

C = BLU - BLU * (E1 / S) * E2LU;

end

