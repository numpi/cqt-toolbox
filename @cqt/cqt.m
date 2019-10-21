classdef cqt
    %CQT  class of continuosly quasi-Toeplitz matrices
    %
    %     T = CQT(neg, pos, A) creates the semi-infinite CQT-matrix with
    %     the specified symbol and finite correction A
    %
    %     T = CQT(neg, pos, U, V) creates the semi-infinite CQT-matrix with
    %     the specified symbol and finite correction U * V.'
    %
    %     T = CQT(neg, pos, A, B, m, n) creates the (m x n)-CQT-matrix with
    %     the specified symbol, top-left correction A and bottom-right 
    %     correction B. If one between m and n is 'inf' then B is ignored.
    %
    %     T = CQT(neg, pos, U, V, W, Z, m, n) creates the 
    %     (m x n)-CQT-matrix with the specified symbol and finite 
    %     corrections A = U * V.' , B = W * Z.'. If one between m and n is
    %     'inf' then B is ignored.
    %
    %     T = CQT(neg, pos) creates the semi-infinite CQT-matrix with the 
    %     specified symbol and an empty finite correction.
    %
    %     T = CQT(A) creates the semi-infinite CQT-matrix with a symbol 
    %     equal to 0  and finite correction equal to A
	%
	%     All the constructors can be called using the syntax:
	%
	%       T = CQT('extended', arg1, ..., argk, v); 
	%
	%     This is equivalent to the commands:
	%
	%       T = CQT(arg1, ..., argk);
	%       T = EXTEND(T, v);
	%
	%     which constructs a quasi-Toeplitz matrix plus a rank 1 correction
	%     of the form ev^T, where e is the vector of all ones. 
	%
	%     T = CQT('hankel', f) constructs the Hankel matrix with symbol
	%     given by f. 
	%
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
        
        % Rank 1 correction containing the limits of the columns
        c
    end
    
    methods
        
        function obj = cqt(varargin)
            %CQT Create a a QuasiToeplitz matrix in the Wiener class
			
			if ~isempty(varargin) && ischar(varargin{1})
				switch varargin{1}
					case 'hankel'
						if length(varargin) < 2
							error('Please specify a symbol for ' + ...
								  'the Hankel matrix');
						end
						
						[U, V] = hankel_compress(varargin{2}, [], ...
							cqtoption('compression'));
						
						obj = cqt([], [], U, V);
						
					case 'extended'
						obj = extend(cqt(varargin{2:end-1}), ...
									 varargin{end});
						
					case 'handle'
						if length(varargin) < 2
							error('Please specify a function to ' + ...
								'evaluate the symbol');
						end
						
						[fm, fp] = evinterp(varargin{2}, ...
							cqtoption('threshold'), [0 1], 0);
						
						obj = cqt(fm, fp);
						
					otherwise
						error('Unsupported constructor option');
				end
				
				return
			end
			
            switch length(varargin)
                case 1
                    obj.n = [];
                    obj.p = [];
                    
                    if size(varargin{1}, 1) <= size(varargin{1}, 2)
                        obj.U = eye(size(varargin{1}, 1));
                        obj.V = varargin{1}';
                    else                    
                        obj.U = varargin{1};
                        obj.V = eye(size(varargin{1}, 2));
                    end
                    
                    obj.sz = [inf,inf];
                    obj.c = zeros(1, 0);
                case 2
                    if ~isempty(varargin{1}) && ...
                            (varargin{1}(1) ~= varargin{2}(1))
                        error([ 'The coefficients of degree 0 do' ...
                            ' not coincide']);
                    end
                    obj.n = varargin{1};
                    obj.p = varargin{2};
                    obj.sz = [inf,inf];
                    obj.c = zeros(1, 0);
                case 3
                    if ~isempty(varargin{1}) && ...
                            (varargin{1}(1) ~= varargin{2}(1))
                        error(['The coefficients of degree 0' ...
                              ' does not coincide' ]);
                    end
                    obj.n = varargin{1};
                    obj.p = varargin{2};
                    
                    if size(varargin{3}, 1) <= size(varargin{3}, 2)
                        obj.U = eye(size(varargin{3}, 1));
                        obj.V = varargin{3}';
                    else                    
                        obj.U = varargin{3};
                        obj.V = eye(size(varargin{3}, 2));
                    end
                    
                    obj.sz = [inf,inf];
                    obj.c = zeros(1, 0);
                case 4
                    if ~isempty(varargin{1}) && ...
                            (varargin{1}(1) ~= varargin{2}(1))
                        error([ 'The coefficients of degree 0' ... 
                              ' does not coincide' ]);
                    end
                    if size(varargin{3},2) == size(varargin{4},2)
                        obj.n = varargin{1};
                        obj.p = varargin{2};
                        obj.U = varargin{3};
                        obj.V = varargin{4};
                        obj.sz = [inf,inf];
                        obj.c = zeros(1, 0);
                    else
                        error(['Incompatible dimensions: the last two' ...
                              ' arguments must have the same number'   ...
                              'of columns' ]);
                    end
                case 6
                    
                    if ~isempty(varargin{1}) && ...
                            (varargin{1}(1) ~= varargin{2}(1))
                        error([ 'The coefficients of degree 0' ...
                                ' does not coincide' ]);
                    end
                    if length(varargin{1}) > varargin{5} || ...
                            length(varargin{2}) > varargin{6}
                        error([ 'Size of the symbol bigger than ' ...
                                'the size of the matrix' ]);
                    end
                    if max(size(varargin{3},1),size(varargin{4},1)) > ...
                            varargin{5} || max(size(varargin{3},2), ...
                            size(varargin{4},2)) > varargin{6}
                        error(['Size of the corrections bigger than'  ...
                              ' the size of the matrix' ]);
                    end
                    if  min([varargin{5},varargin{6}]) < 0
                        error(['Invalid parameter for the size' ...
                              ' of the corrections' ]);
                    end
                    obj.n = varargin{1};
                    obj.p = varargin{2};
                    
                    if size(varargin{3}, 1) <= size(varargin{3}, 2)
                        obj.U = eye(size(varargin{3}, 1));
                        obj.V = varargin{3};
                    else                    
                        obj.U = varargin{3};
                        obj.V = eye(size(varargin{3},2));
                    end
                    
                    if max (varargin{5},varargin{6}) ~=inf
                        if size(varargin{4}, 1) <= size(varargin{4}, 2)
                            obj.W = eye(size(varargin{4}, 1));
                            obj.Z = varargin{4}(end:-1:1, end:-1:1);
                        else
                            obj.W = varargin{4}(end:-1:1,end:-1:1);
                            obj.Z = eye(size(varargin{4},2));
                        end
                    elseif ~isempty(varargin{4})
                        warning(['The bottom right correction is' ...
                            'ignored due to the infinite dimension ' ...
                            'of the matrix' ]);
                    end
                    obj.sz = [varargin{5},varargin{6}];
                    obj.c = zeros(1, 0);
                    obj = merge_corrections(obj);
                case 8
                    if ~isempty(varargin{1}) && ...
                            (varargin{1}(1) ~= varargin{2}(1))
                        error(['The coefficients of degree 0' ...
                            ' do not coincide' ]);
                    end
                    if max(size(varargin{3},1),size(varargin{5},1)) > ...
                            varargin{7} || ...
                            max(size(varargin{4},1), ... 
                            size(varargin{6},1)) > varargin{8}
                        error(['Size of the corrections bigger than' ... 
                            ' the size of the matrix']);
                    end
                    if min([varargin{7},varargin{8}]) < 0
                        error(['Invalid parameter for the size'  ...
                            ' of the corrections']);
                    end
                    if size(varargin{3},2) ~= size(varargin{4},2) || ...
                            size(varargin{5},2) ~= size(varargin{6},2)
                        error([ 'Factors of the corrections must have' ...
                            'the same number of columns' ]);
                    end
                    obj.n = varargin{1};
                    obj.p = varargin{2};
                    obj.U = varargin{3};
                    obj.V = varargin{4};
                    if max (varargin{7},varargin{8}) ~= inf
                        obj.W = varargin{5}(end:-1:1,end:-1:1);
                        obj.Z = varargin{6}(end:-1:1,end:-1:1);
                    elseif max(length(varargin{5}),length(varargin{6})) > 0
                        warning([ 'The bottom right correction is'   ...
                            ' ignored due to the infinite dimension' ...
                            ' of the matrix' ]);
                    end
                    obj.sz = [varargin{7},varargin{8}];
                    obj.c = zeros(1, 0);
                    obj = merge_corrections(obj);
                otherwise
                    error('Invalid number of parameters');
            end
            
            if (~isempty(obj.n) && ~isvector(obj.n)) || ...
                    (~isempty(obj.p) && ~isvector(obj.p))
                error('The symbol must be a 1D vector');
            end
            
            % Make sure that the symbol is represented as a row vector.
            obj.n = reshape(obj.n, 1, length(obj.n));
            obj.p = reshape(obj.p, 1, length(obj.p));
        end
        
        function A = merge_corrections(A)
            %MERGE_CORRECTIONS
            m = size(A, 1);
            n = size(A, 2);
            
            if size(A.U, 1) + size(A.W, 1) > n && ...
                    size(A.V, 1) + size(A.Z, 1) > n
                A.U = [ [ A.U ; zeros(m-size(A.U,1),size(A.U,2)) ] , ...
                    [ zeros(m-size(A.W,1),size(A.W,2)) ; ...
                    A.W(end:-1:1,end:-1:1) ] ];
                A.V = [ [ A.V ; zeros(n-size(A.V,1),size(A.V,2)) ] , ...
                    [ zeros(n-size(A.Z,1),size(A.Z,2)) ; ...
                    A.Z(end:-1:1,end:-1:1) ] ];
                A.W = [];
                A.Z = [];
            end
        end
        
    end
end
