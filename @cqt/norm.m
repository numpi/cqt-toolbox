function r = norm(T, p)
%NORM CQT-norm of a CQT-matrix
%
%     r = NORM(T) computes the CQT-norm of a CQT-matrix

if max(T.sz) == inf
	if exist('p', 'var') && ~(strcmp(p, 'CQT') || strcmp(p, 'cqt'))
		error('Only the CQT norm is supported for infinite matrices');
	end
	
	r = qt_norm(T.n, T.p, T.U, T.V);
else
	if max(T.sz) == inf
		T.sz = [ min(max(size(T.U,1), length(T.n) + T.sz(2)), T.sz(1)), ...
			min(max(size(T.V,2), length(T.p) + T.sz(1)), T.sz(2)) ];
	end
	% For small matrices, or the ones that have overlapping corrections, we
	% compute the norm of the full version
	m = size(T, 1);
	n = size(T, 2);
	if max(T.sz) < 50 || ((size(T.U,1) + size(T.W,1) >= m) && (size(T.V,1) + size(T.Z,1) >= n) )
		if exist('p', 'var')
			r = norm(full(T), p);
		else
			r = norm(full(T));
		end
	else
		if ~exist('p', 'var')
			p = 2;
		end
		
		switch p
			case 1
				r = fqt_norm_1(T);
			case 2
				r = fqt_norm_2(T);
			case inf
				r = fqt_norm_inf(T);
			case 'cqt'
				r = qt_norm(T.n, T.p, T.U, T.V) + qt_norm(0, 0, T.W, T.Z);
			otherwise
				error('Unsupported norm');
		end
	end
	
	
end

