function varargout = size(T, idx)
%SIZE Obtain the size of the CQT matrix.
if nargout == 1
	sz = T.sz;
	
	if exist('idx', 'var')
		if idx < 1 || idx > 2
			error('Invalid dimension specified');
		else
			sz = sz(idx);
		end
	end
	varargout{1} = sz;
else
	varargout{1} = T.sz(1);
	varargout{2} = T.sz(2);
end

