function X = mldivide(A, B)

if isa(B, 'cqt') && ~isa(A, 'cqt')
	if isscalar(A)
		X = (1 / A) * B;
		return;
	end
	
	if B.sz(1) ~= size(A,2)
		error('Incompatible dimensions');
	end
	cl = B.sz(2);
	if cl == inf
		B.sz(2) = max(size(B.V,1), length(B.p) + B.sz(1));
	end
	
	X = cqt([], [], A \ [B.U;zeros(size(A,1)-size(B.U,1),size(B.U,2))], B.V,...
		A \ [zeros(size(A,1)-size(B.W,1),size(B.W,2));B.W(end:-1:1,end:-1:1)], B.Z(end:-1:1,end:-1:1), size(A,1), cl);
	if ~isempty(B.n) && ~isempty(B.p)
		X = X + cqt([],[], A\toep([B.n],[B.p],size(A,2), length(B.p)+size(A,2)-1),[], X.sz(1),X.sz(2));
	end
	return;
end

if isa(A, 'cqt') && ~isa(B, 'cqt')
	if A.sz(2) ~= size(B, 1)
		error('Incompatible dimensions');
	end
	
	cB = cqt([], [], B, [], size(B, 1), size(B, 2));
	X  = full(A \ cB);
	return;
end

if isa(A,'cqt') && isa(B,'cqt')
	if isinf(A.sz(1))
		X = qt_mldivide(A, B);
	else
		X = fqt_mldivide(A, B);
	end
end
