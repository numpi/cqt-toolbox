function I = eye(varargin)
%EYE Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3 && ischar(varargin{2}) && strcmp(varargin{2}, 'like')
    if ~isa(varargin{3}, 'cqt')
        error('Unsupported call to cqt/eye');
    end
    
    n = varargin{1}(1);
    
    I = cqt(1, 1, [], [], [], [], n, n);
end

end

