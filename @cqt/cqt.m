classdef cqt 
%CQT  class of continuosly quasi-Toeplitz matrices
%
%     T = CQT(pos, neg, A) creates the semi-infinite CQT-matrix with 
%     the specified symbol and finite correction A
%
%     T = CQT(pos, neg, U, V) creates the semi-infinite CQT-matrix with 
%     the specified symbol and finite correction U * V.'
%
%     T = CQT(pos, neg, A, B, m, n) creates the finite CQT-matrix with 
%     the specified symbol and finite corrections A , B
%
%     T = CQT(pos, neg, U, V, W, Z, m, n) creates the finite CQT-matrix with 
%     the specified symbol and finite corrections U * V.' , W * Z.'
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
				if varargin{1}(1) ~= varargin{2}(1)
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
	    			if varargin{1}(1) ~= varargin{2}(1)
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
				if varargin{1}(1) ~= varargin{2}(1)
					error('The coefficients of degree 0 does not coincide');
	    			end
				if length(varargin{1}) > varargin{5} || length(varargin{2}) > varargin{6}
					error('Size of the symbol bigger than the size of the matrix');
	    			end
				if max(size(varargin{3},1),size(varargin{4},1)) > varargin{5} ||...
			 		max(size(varargin{3},2),size(varargin{4},2)) > varargin{6}
					error('Size of the corrections bigger than the size of the matrix');
				end
				if max([varargin{3},varargin{4}]) == inf || min([varargin{3},varargin{4}]) <= 0
					error('Invalid parameter for the size of the corrections');
				end
				obj.n = varargin{1};
            			obj.p = varargin{2};
            			obj.U = varargin{3};
            			obj.V = eye(size(varargin{3},2));
				obj.W = varargin{4}(end:-1:1,end:-1:1);
				obj.Z = eye(size(varargin{4},2));
				obj.sz = [varargin{5},varargin{6}];
			case 8
				if varargin{1}(1) ~= varargin{2}(1)
					error('The coefficients of degree 0 do not coincide');
	    			end
				if max(size(varargin{3},1),size(varargin{5},1)) > varargin{7} ||...
			 		max(size(varargin{4},1),size(varargin{6},1)) > varargin{8}
					error('Size of the corrections bigger than the size of the matrix');
				end
				if max([varargin{7},varargin{8}]) == inf || min([varargin{7},varargin{8}]) <= 0
					error('Invalid parameter for the size of the corrections');
				end
				if size(varargin{3},2) ~= size(varargin{4},2) || size(varargin{5},2) ~= size(varargin{6},2)
					error('Factors of the corrections must have the same number of columns');
				end
				obj.n = varargin{1};
            			obj.p = varargin{2};
            			obj.U = varargin{3};
            			obj.V = varargin{4};
				obj.W = varargin{5}(end:-1:1,end:-1:1);
				obj.Z = varargin{6}(end:-1:1,end:-1:1);
				obj.sz = [varargin{7},varargin{8}];
			otherwise
				error('Invalid number of parameters');
        	end
    	end
     end
end
