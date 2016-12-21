function [varargout] = correction(T)
%CORRECTION Returns the finite correction to the semiinfinite Toeplitz T. 
%
% E = CORRECTION(T) returns the finite correction to the Toeplitz part of
%     the CQT matrix T. 
%
% [U, V] = CORRECTION(T) returns a low-rank factorization of the matrix E
%     obtained calling E = CORRECTION(T), that is E = U * V.'. 

if nargout <= 1
    varargout{1} = T.U * T.V.';
elseif nargout == 2
    varargout{1} = T.U;
    varargout{2} = T.V;
else
    error('Unsupported number of output arguments');
end

end

