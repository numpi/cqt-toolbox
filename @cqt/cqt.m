classdef cqt 
%CQT  class of continuosly quasi-Toeplitz matrices
%
%     T = CQT(pos, neg, U, V) creates the CQT-matrix with the specified symbol
%     and finite correction U * V.'
%
%     T = CQT(pos, neg) creates the CQT-matrix with the specified symbol and
%     an empty finite correction
%
%     T = CQT(A) creates the CQT-matrix with a symbol equal to 0 and finite
%     correction equal to A
    properties
        % Vector containing the nonnegative coefficients of the symbol
        p
        
        % Vector containing the nonpositive coefficients of the symbol
        n
        
        % Left factor of the finite top-left correction
        U
        
        % Right factor of the finite top-left correction
        V
    end
    
    methods

        function obj = cqt(varargin)
        %CQT Create a a QuasiToeplitz matrix in the Wiener class
		if length(varargin) == 4
	    		if varargin{1}(1) ~= varargin{2}(1)
				error('The coefficients of degree 0 does not coincide');
	    		end
	    		if size(varargin{3},2) == size(varargin{4},2)
            			obj.n = varargin{1};
            			obj.p = varargin{2};
            			obj.U = varargin{3};
            			obj.V = varargin{4};
	    		else
				error('Incompatible dimensions: the last two arguments must have the same number of columns');
	   		end
		elseif length(varargin) == 1
	    		obj.n = 0;
            		obj.p = 0;
            		obj.U = varargin{1};
            		obj.V = eye(size(varargin{1},2));
		elseif length(varargin) == 2
			if varargin{1}(1) ~= varargin{2}(1)
				error('The coefficients of degree 0 does not coincide');
	    		end
	    		obj.n = varargin{1};
            		obj.p = varargin{2};
        	end

    	end
     end
end
