function value = getoption(keyword)
%GETOPTION Obtain the value for a global option

global cqt_inversion;
global cqt_sqrt;
global cqt_compression;
global cqt_threshold;
global cqt_wiener_hopf;

switch keyword
    case 'inversion'
        value = cqt_inversion;
        
    case 'sqrt'
        value = cqt_sqrt;

    case 'wiener-hopf'
        value = cqt_wiener_hopf;
        
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

