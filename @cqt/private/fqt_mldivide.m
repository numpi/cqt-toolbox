function C = fqt_mldivide(A, B)
%FQT_MLDIVIDE

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
	
	if isa(B, 'cqt') || max(length(A.n), length(A.p)) > 8
		C = inv(A) * B;
	else
		C = fqt_solvetriang(A, B);
	end
	
	return;
end

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, E.W(end:-1:1,end:-1:1), A.sz(1), size(E.U,2) + size(E.W,2));
E2 = cqt([], [], E.V, E.Z(end:-1:1,end:-1:1), A.sz(2), size(E.V,2) + size(E.Z,2));

% Compute A \ B by Sherman-Morrison
Linv = inv(L);
Uinv = inv(U);
LUB = Linv * (Uinv * B);
LUE1 = Linv * (Uinv * E1);

S = eye(size(E1, 2)) + full(E2.' * LUE1);

C = LUB - LUE1 * (S \ (E2.' * LUB));