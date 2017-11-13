function A = extend(A, c)
%EXTEND Add a rank-1 e * c' correction to A, where e = ones(inf,1).
%
% If the correction has already been set, the new value of C replaces it.

if min(size(A)) ~= inf
    error('extended Toeplitz arithmetic is only supported for infinite matrices');
end

A.c = c;

end

