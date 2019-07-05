function varargout = cqtgallery(name, varargin)
%CQTGALLERY Construct a test example, similar to MATLAB's gallery.
%
% [out1, out2, ...] = CQTGALLERY(name, param1, param2, ...) takes a name 
% for a given example and generates some CQT matrices representing the 
% given problem. The meaning of the parameters and the output is
% problem-dependent.
%
% Available problems:
%
%  jackson     QBD processes arising in waiting queues [1], taken from [2].
%              There are 10 different test problems available, which can
%              be obtained by calling CQTGALLERY('jackson', N), where N is
%              in [1, ..., 10].
%
% [1] J. R. Jackson. Networks of waiting lines. Operations research,
%     5(4):518–521, 1957.
%
% [2] A. J. Motyer and P. G. Taylor. Decay rates for quasi-birth-and-death
%     processes with countablu many phases and tridiagonal block
%     generators. Adv. Appl. Prob., 38:522–544, 2006.

switch name
    case 'jackson'
        if length(varargin) < 1 || ~isnumeric(varargin{1})
            error(...
                'The argument for the Jackson problem has to be numeric');
        end
        [Am1, A0, A1, hA0] = jackson(varargin{1});
        
        varargout = {Am1, A0, A1, hA0};
	case 'rand'
		if length(varargin) >= 1 && ~isnumeric(varargin{1})
			error('Invalid argument to the rand constructor');
		end
		
		am = rand(1, 5);
		ap = rand(1, 5); ap(1) = am(1);
		
		varargout{1} = cqt(am, ap, rand(10));
    otherwise
        error('Invalid problem name specified');
end

end

