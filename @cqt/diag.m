function D = diag(C)
%DIAG Extract the diagonal of a CQT matrix.
%
% D = DIAG(C) obtains a vector with the diagonal of the CQT matrix C.
% Notice this only makes sense for finite matrices.

if min(size(C)) == inf
    error('DIAG can only be used on finite CQT matrices');
end

D = C;

if ~isempty(D.p)
    D.p = D.p(1);
    D.n = D.n(1);
end

D.U = [];
D.V = [];
D.W = [];
D.Z = [];

t = C.U * C.V.';

% This workaround is needed for the case where t is 1 x N or N x 1,
% otherwise diag creates a diagonal matrix instead of extracting the
% diagonal entries.
if min(size(t)) == 1
    td = t(1);
else
    td = diag(C.U * C.V.');
end

D.U = diag(td);
D.V = eye(size(D.U));

if max(C.sz) < inf
    t = C.W * C.Z.';
    
    % Same trick as above
    if min(size(t)) == 1
        td = t(1);
    else
        td = diag(C.U * C.V.');
    end
    
    D.W = diag(td);
    D.Z = eye(size(D.W));
end


end

