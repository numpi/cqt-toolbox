function setoption(keyword, value)
%SETOPTION Set an option for the CQT toolbox.

global cqt_inversion;
global cqt_sqrt;
global cqt_compression;
global cqt_threshold;
global cqt_wiener_hopf;

switch keyword
    case 'inversion'
        cqt_inversion = value;
    case 'sqrt'
        cqt_sqrt = value;
    case 'wiener-hopf'
        cqt_wiener_hopf = value;
    case 'compression'
        cqt_compression = value;
    case 'threshold'
        cqt_threshold = value;
    otherwise
        error('Unsupported option specified');
end

end

