function value = getoption(keyword)
%GETOPTION Obtain the value for a global option

global cqt_inversion;
global cqt_sqrt;

switch keyword
	case 'inversion'
		value = cqt_inversion;
		
	case 'sqrt'
		value = cqt_sqrt;
		
	otherwise
		error('Unsupported option specified');
end

