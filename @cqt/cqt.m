classdef cqt
	%CQT  class of continuosly quasi-Toeplitz matrices
	%
	%     T = CQT(pos, neg, A) creates the semi-infinite CQT-matrix with
	%     the specified symbol and finite correction A
	%
	%     T = CQT(pos, neg, U, V) creates the semi-infinite CQT-matrix with
	%     the specified symbol and finite correction U * V.'
	%
	%     T = CQT(pos, neg, A, B, m, n) creates the (m x n)-CQT-matrix with
	%     the specified symbol, top-left correction A and bottom-right correction B.
	%     If one between m and n is 'inf' then B is ignored.
	%
	%     T = CQT(pos, neg, U, V, W, Z, m, n) creates the (m x n)-CQT-matrix with
	%     the specified symbol and finite corrections A = U * V.' , B = W * Z.'.
	%     If one between m and n is 'inf' then B is ignored.
	%
	%     T = CQT(pos, neg) creates the semi-infinite CQT-matrix with the specified
	%     symbol and an empty finite correction
	%
	%     T = CQT(A) creates the semi-infinite CQT-matrix with a symbol equal to 0
	%     and finite correction equal to A
	properties
		% Vector containing the nonnegative coefficients of the symbol
		p
		
		% Vector containing the nonpositive coefficients of the symbol
		n
		
		% Left factor of the finite top-left correction
		U
		
		% Right factor of the finite top-left correction
		V
		
		% Left factor of the finite bottom-left correction
		W
		
		% Right factor of the finite right-left correction
		Z
		
		% Vector containing the size of the CQT-matrix
		sz
	end
	
	methods
		
		function obj = cqt(varargin)
			%CQT Create a a QuasiToeplitz matrix in the Wiener class
			switch length(varargin)
				case 1
					obj.n = 0;
					obj.p = 0;
					obj.U = varargin{1};
					obj.V = eye(size(varargin{1},2));
					obj.sz = [inf,inf];
				case 2
					if length(varargin{1}) > 0 && (varargin{1}(1) ~= varargin{2}(1))
						error('The coefficients of degree 0 does not coincide');
					end
					obj.n = varargin{1};
					obj.p = varargin{2};
					obj.sz = [inf,inf];
				case 3
					if varargin{1}(1) ~= varargin{2}(1)
						error('The coefficients of degree 0 does not coincide');
					end
					obj.n = varargin{1};
					obj.p = varargin{2};
					obj.U = varargin{3};
					obj.V = eye(size(varargin{3},2));
					obj.sz = [inf,inf];
				case 4
					if length(varargin{1}) > 0 && (varargin{1}(1) ~= varargin{2}(1))
						error('The coefficients of degree 0 does not coincide');
					end
					if size(varargin{3},2) == size(varargin{4},2)
						obj.n = varargin{1};
						obj.p = varargin{2};
						obj.U = varargin{3};
						obj.V = varargin{4};
						obj.sz = [inf,inf];
					else
						error('Incompatible dimensions: the last two arguments must have the same number of columns');
					end
				case 6
					
					if length(varargin{1}) > 0 && (varargin{1}(1) ~= varargin{2}(1))
						error('The coefficients of degree 0 does not coincide');
					end
					if length(varargin{1}) > varargin{5} || length(varargin{2}) > varargin{6}
						error('Size of the symbol bigger than the size of the matrix');
					end
					if max(size(varargin{3},1),size(varargin{4},1)) > varargin{5} ||...
							max(size(varargin{3},2),size(varargin{4},2)) > varargin{6}
						error('Size of the corrections bigger than the size of the matrix');
					end
					if  min([varargin{5},varargin{6}]) < 0
						error('Invalid parameter for the size of the corrections');
					end
					obj.n = varargin{1};
					obj.p = varargin{2};
					obj.U = varargin{3};
					obj.V = eye(size(varargin{3},2));
					if max (varargin{5},varargin{6}) ~=inf
						obj.W = varargin{4}(end:-1:1,end:-1:1);
						obj.Z = eye(size(varargin{4},2));
					elseif length(varargin{4}) > 0
						warning('The bottom right correction is ignored due to the infinite dimension of the matrix');
					end
					obj.sz = [varargin{5},varargin{6}];
				case 8
					if length(varargin{1}) > 0 && (varargin{1}(1) ~= varargin{2}(1))
						error('The coefficients of degree 0 do not coincide');
					end
					if max(size(varargin{3},1),size(varargin{5},1)) > varargin{7} ||...
							max(size(varargin{4},1),size(varargin{6},1)) > varargin{8}
						error('Size of the corrections bigger than the size of the matrix');
					end
					if min([varargin{7},varargin{8}]) < 0
						error('Invalid parameter for the size of the corrections');
					end
					if size(varargin{3},2) ~= size(varargin{4},2) || size(varargin{5},2) ~= size(varargin{6},2)
						error('Factors of the corrections must have the same number of columns');
					end
					obj.n = varargin{1};
					obj.p = varargin{2};
					obj.U = varargin{3};
					obj.V = varargin{4};
					if max (varargin{7},varargin{8}) ~= inf
						obj.W = varargin{5}(end:-1:1,end:-1:1);
						obj.Z = varargin{6}(end:-1:1,end:-1:1);
					elseif max(length(varargin{5}),length(varargin{6})) > 0
						warning('The bottom right correction is ignored due to the infinite dimension of the matrix');
					end
					obj.sz = [varargin{7},varargin{8}];
				otherwise
					error('Invalid number of parameters');
			end
			
			if (~isempty(obj.n) && ~isvector(obj.n)) || (~isempty(obj.p) && ~isvector(obj.p))
				error('The symbol must be a 1D vector');
			end
			
			% Make sure that the symbol is represented as a row vector.
			obj.n = reshape(obj.n, 1, length(obj.n));
			obj.p = reshape(obj.p, 1, length(obj.p));
		end
	end
end
