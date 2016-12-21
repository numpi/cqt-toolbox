classdef cqt 
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
            obj.n = varargin{1};
            obj.p = varargin{2};
            obj.U = varargin{3};
            obj.V = varargin{4};
        end

    end
end