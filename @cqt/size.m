function varargout = size(T, idx)
%SIZE Obtain the size of the CQT matrix.

if exist('idx', 'var')
    if idx < 1 || idx > 2
        error('Invalid dimension specified');
    else
        sz = T.sz(idx);
    end
    varargout{1} = sz;
else
    if nargout == 2
        varargout{1} = T.sz(1);
        varargout{2} = T.sz(2);
    else
        varargout{1} = [ T.sz(1), T.sz(2) ];
    end
end

