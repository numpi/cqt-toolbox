function M = slicemat_consecutive(T, imin, imax, jmin, jmax)
%SLICEMAT_CONSECUTIVE Slice a contiguous part of a CQT matrix.
%
%     M = SLICEMAT_CONSECUTIVE(T, IMIN, IMAX, JMIN, JMAX) extracts a dense
%     representation of a slice of a CQT matrix with row indices between
%     IMIN and IMAX and column indices between JMIN and JMAX.

M = zeros(imax - imin + 1, jmax - jmin + 1);

% Add the correction at the top-left corner
if imin <= size(T.U, 1) && jmin <= size(T.V, 1)
	max_row_index = min(imax, size(T.U, 1));
	max_col_index = min(jmax, size(T.V, 1));
	
	max_rel_row_index = max_row_index - imin + 1;
	max_rel_col_index = max_col_index - jmin + 1;
	
	M(1:max_rel_row_index,1:max_rel_col_index) = ...
		M(1:max_rel_row_index,1:max_rel_col_index) + ...
		T.U(imin:max_row_index,:) * ...
		T.V(jmin:max_col_index,:).';
end

% Add the correction at the bottom-right corner
if imax >= T.sz(1) - size(T.W,1) && jmax >= T.sz(2) - size(T.Z,1)
	min_row_index = max(imin, T.sz(1) - size(T.W, 1) + 1);
	min_col_index = max(jmin, T.sz(2) - size(T.Z, 1) + 1);
	
	min_rel_row_index = min_row_index - imin + 1;
	min_rel_col_index = min_col_index - jmin + 1;	
	
	row_idx = (T.sz(1) - min_row_index + 1) : -1 : (T.sz(1) - imax + 1);
	col_idx = (T.sz(2) - min_col_index + 1) : -1 : (T.sz(2) - jmax + 1);
	
	M(min_rel_row_index:end,min_rel_col_index:end) = ...
		M(min_rel_row_index:end,min_rel_col_index:end) + ...
		T.W(row_idx, :) * T.Z(col_idx, :).';
end

% Add the symbol
if ~isempty(T.n)
	mp = length(T.n);
	smb = [ T.n(end:-1:1), T.p(2:end) ];
	
	% Compute the new middle point
	mp = mp + jmin - imin;

	if mp > jmin - jmax - 1 && mp < length(smb) + imax - imin + 1
		if mp <= 0
			neg = 0;
			pos = [ zeros(1, 1 - mp) , smb ];
		end
		
		if mp > 0 && mp <= length(smb)
			neg = smb(mp:-1:1);
			pos = smb(mp:end);			
		end
		
		if mp > length(smb)
			pos = 0;
			neg = [ zeros(1, mp - length(smb)), smb(end:-1:1) ];
		end
		
		M = M + toep(neg, pos, size(M, 1), size(M, 2));
	end
end