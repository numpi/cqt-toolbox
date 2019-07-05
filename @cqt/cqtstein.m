function X = cqtstein(A, B, C, varargin)
%CQTSTEIN Solves the equation A*X*B + X = C
%
% X = CQTSTEIN(A, B, C) computes the CQT solution to the matrix equation 
%
%    (1)   A*X*B + X + C = 0
%
%     under the assumption that both A and B have spectral radius strictly
%     smaller than 1. 

p = inputParser;

addParameter(p, 'tol', cqtoption('threshold'));
addParameter(p, 'method', 'galerkin');
addParameter(p, 'debug', false);

parse(p, varargin{:});

tol = p.Results.tol;
method = p.Results.method;
debug = p.Results.debug;

maxit = 1000;

if strcmp(method, 'fixedpoint')
	X = cqt([], []);
	for j = 1 : maxit
		Xold = X;
		X = -A*X*B - C;
		
		res = norm(X - Xold, inf) / norm(X, inf);
		if debug
			fprintf('FIXEDPOINT :: Iteration %d, res = %e\n', j, res);
		end
		
		if res < tol
			break;
		end
	end
	return;
end

IA = eye(size(A), 'like', A);
IB = eye(size(B), 'like', B);

% M1 = inv(A - IA);
% M2 = inv(B + IB);

M1 = inv(IA - A);
M2 = inv(B + IB);

switch method
    case 'galerkin'
        X = cqtlyap(M1 * (A + IA), M2 * (IB - B), 2 * M1 * C * M2, ...
            'tol', tol, 'debug', debug, 'poles', [-ones(1,11) inf ]);
    case 'adi'
        s = max(norm(A, inf), norm(B, inf));

        b = (1 - s) / (1 + s);
		% b = 1;
        
        % Alternative choice for the parameters that does not depend on the
        % matrix is a = -1, b = 1. This always works, but convergence is
        % slower. 
        
		b = -1;
        a = -b;
        
		X = cqtadi(M1 * (A + IA), M2 * (IB - B), 2 * M1 * C * M2, a, b, ...
            'tol', tol, 'debug', debug);
		%norm(M1 * (A + IA) * X + X * M2 * (IB - B) + 2 * M1 * C * M2)
		%norm(X)
    otherwise
        error('Unsupported Sylvester solver selected');
end

end