function value = getoption(keyword)
%GETOPTION Obtain the value for a global option

global inversion;

switch keyword     
    case 'inversion'
        value = inversion;
        
    otherwise
        error('Unsupported option specified');
end

