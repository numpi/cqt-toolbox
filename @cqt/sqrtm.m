function X = sqrtm(A)
%SQRTM Compute the matrix square root of X.
%
% Y = SQRT(X) computes a matrix Y such that Y^2 = X. The computed square
% root is the principal branch of the square root.

converged = false;
threshold = 1.0e2 * eps * norm(A, inf);

switch cqtoption('sqrt')
	
	case 'db'
		
		X = A;
		Y = cqt(1, 1, [], [], A.sz(1), A.sz(2));
		
		while ~converged
			Xnew = .5 * (X + inv(Y));
			Ynew = .5 * (Y + inv(X));
			
			if norm(Xnew - X, inf) < threshold
				converged = true;
			end
			
			X = Xnew;
			Y = Ynew;
		end
		
	case 'cr'
		X = A;
		E = .5 * (cqt(1, 1, [], [], A.sz(1), A.sz(2)) - A);
		
		while ~converged
			Xnew = X + E;
			Enew = - .5 * (E * inv(Xnew) * E);
			
			if norm(Enew, inf) < threshold
				converged = true;
			end
			
			X = Xnew;
			E = Enew;
		end
		
	otherwise
		error('Unsupported value for the option "sqrt"');
end


end

