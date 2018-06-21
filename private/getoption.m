function value = getoption(keyword)
%GETOPTION Obtain the value for a global option

global cqt_inversion;
global cqt_sqrt;
global cqt_compression;
global cqt_threshold;

switch keyword
    case 'inversion'
        value = cqt_inversion;
        
    case 'sqrt'
        value = cqt_sqrt;
        
    case 'compression'
        if isempty(cqt_compression)
            cqt_compression = 'lanczos';
        end
        value = cqt_compression;
        
    case 'threshold'
        if isempty(cqt_threshold)
            cqt_threshold = 1e-12;
        end
        value = cqt_threshold;
        
    otherwise
        error('Unsupported option specified');
end

