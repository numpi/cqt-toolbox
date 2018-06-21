function varargout = cqtrank(C)
%CQTPART Rank of the correction to the Toeplitz part.
%
% K = CQTRANK(C) obtains the rank of the correction to the Toeplitz part in
% C. This is done account for both top and bottom corrections for finite
% matrices.
%
% [KT, KB] = CQTRANK(C) retrieves the ranks of the top and bottom
% corrections separately.

if nargout == 1
    varargout{1} = size(C.U, 2) + size(C.W, 2);
else
    varargout{1} = size(C.U, 2);
    varargout{2} = size(C.W, 2);
end

end

