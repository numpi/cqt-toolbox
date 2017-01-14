function sz = size(T, idx)
%SIZE Obtain the size of the CQT matrix.

sz = T.sz;

if exist('idx', 'var')
	if idx < 1 || idx > 2
		error('Invalid dimension specified');
	else
		sz = sz(idx);
	end
end

end

