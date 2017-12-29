function optval = cqtoption(keyword, value)
%CQTOPTION Set or get an option for the CQT toolbox.
%
% CQTOPTION(KEYWORD, VALUE) sets the option KEYWORD to the specified value.
% The possible options are the following:
%
%   'inversion': [ 'cr', 'fft' ]
%      Select the algorithm used to perform the Toeplitz inversion.
%
%   'sqrt': [ 'db', 'cr' ]
%      Select the iteration used to compute the square root.
%
%   'threshold': 
%      The threshold used in truncation.
%

switch keyword
	case 'inversion'
		if nargin == 1
			optval = getoption('inversion');
			if isempty(optval)
				optval = 'cr';
			end
		else
			setoption('inversion', value);
		end
		
	case 'sqrt'
		if nargin == 1
			optval = getoption('sqrt');
			if isempty(optval)
				optval = 'db';
			end
		else
			setoption('sqrt', value);
        end
        
    case 'compression'
        if nargin == 1
            optval = getoption('compression');
        else
            setoption('compression', value);
        end
        
    case 'threshold'
        if nargin == 1
            optval = getoption('threshold');
        else
            setoption('threshold', value);
        end
        
end

end

