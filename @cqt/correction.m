function [varargout] = correction(T)
%CORRECTION Returns the finite correction to the semiinfinite Toeplitz T. 
%
% E = CORRECTION(T) returns the finite correction to the Toeplitz part of
%     the CQT matrix T. 
%
% [E, F] = CORRECTION(T) returns the correction E of the top part and the
%     correction F of the trailing block of the matrix. Notice that this is
%     only meaningful if the matrix is finite. 
%
% [U, V] = CORRECTION(T, 'factored', true) returns a low-rank factorization 
%     of the matrix E obtained calling E = CORRECTION(T), that is 
%     E = U * V.'.
%
% [U, V, W, Z] = CORRECTION(T) is similar to the previous
%     form and returns the factored form U * V.' and W * Z.' of the leading
%     and trailing corrections to the Toeplitz matrix. This is only
%     meaningful in the finite case. 
%
% Note: This function has two different behavior depending on the size of
% the matrix. If T is infinite then [U, V] = CORRECTION(T) is the low-rank
% factorization of E, whilst when T is finite the two outputs are assigned
% the dense corrections to the Toeplitz part.

switch T.sz(1)
    case inf
        if nargout <= 1
            varargout{1} = T.U * T.V.';
        else
            varargout{1} = T.U;
            varargout{2} = T.V;
        end
        
    otherwise        
        if nargout <= 2
            varargout{1} = T.U * T.V.';
            varargout{2} = T.W * T.Z.';
            
            % The trailing block is stored in reversed form
            varargout{2} = varargout{2}(end:-1:1, end:-1:1);
        else
            varargout{1} = T.U;
            varargout{2} = T.V;
            varargout{3} = T.W(end:-1:1,end:-1:1);
            varargout{4} = T.Z(end:-1:1,end:-1:1);
        end

end

end