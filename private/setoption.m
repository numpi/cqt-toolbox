function setoption(keyword, value)
%SETOPTION Set an option for the CQT toolbox. 

global inversion;

switch keyword
    case 'inversion'        
        inversion = value;        
    otherwise
        error('Unsupported option specified');
end
        
end

