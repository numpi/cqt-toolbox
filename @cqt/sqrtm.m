function Y = sqrtm(A)
%SQRTM Compute the matrix square root of X. 
%
% Y = SQRT(X) computes a matrix Y such that Y^2 = X. The computed square
% root is the principal branch of the square root. 

Y = A;
converged = false;
threshold = 1.0e2 * eps * norm(A, inf);

while ~converged
	oldY = Y;
	Y = .5 * (Y + Y \ A);
	
	if norm(oldY - Y, inf) < threshold
		converged = true;
	end
end


end

