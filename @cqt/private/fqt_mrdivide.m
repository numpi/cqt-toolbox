function C = fqt_mrdivide(B, A)
%FQT_MRDIVIDE Compute B * inv(A) for finite CQT matrices

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, E.W(end:-1:1,end:-1:1), A.sz(1), size(E.U,2) ...
    + size(E.W,2));
E2 = cqt([], [], E.V, E.Z(end:-1:1,end:-1:1), A.sz(2), size(E.V,2) ...
    + size(E.Z,2));

% Compute A \ B by Sherman-Morrison
Linv = inv(L);
Uinv = inv(U);
BLU = (B * Linv) * Uinv;
E2LU = (E2.' * Linv) * Uinv;

%LU = inv(L)*inv(U);
%BLU = B * LU;
%E2LU = E2.' * LU;

S = eye(size(E1, 2)) + full(E2LU * E1);

SE2LU = cqt([],[], S(1:size(E2LU.U,2),1:size(E2LU.U,2))\E2LU.U, ...
    E2LU.V, S(size(E2LU.U,2)+1:end,size(E2LU.U,2)+1:end) \ ...
    E2LU.W(end:-1:1,end:-1:1), ...
    E2LU.Z(end:-1:1,end:-1:1), E2LU.sz(1), E2LU.sz(2));

E1SE2LU =  -E1 * SE2LU;
E1SE2LU.p = 1;
E1SE2LU.n = 1;
C = BLU * E1SE2LU;

end

