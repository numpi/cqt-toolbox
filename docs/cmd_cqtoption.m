%% Global options for the QT toolbox 
%
% A few options can be set globally in the QT toolbox through the
% |cqtoption| command. 
%
% The supported options are the following:
%
% * |'threshold'|: The truncation threshold in the operations. This is
%   indicative of the accuracy that will be attained by the methods. The
%   default value is |1e-12|.
% * |'inversion'|: Method used for the inversion of the symbol (and
%   therefore of QT matrices); can be either |'cr'| (for Cyclic Reduction) or
%   |'fft'| (Fast Fourier Transform). 
% * |'sqrt'|: Method used to compute the matrix square root. Can be either
%   |'db'| (Denman-Beavers) or |'cr'| (Cyclic Reduction).
% * |'compression'|: The compression algorithm for the compact correction;
%   can be either |'lanczos'| for the Golub-Kahan-Lanczos algorithm for the
%   SVD, or |'random'| to use a randomized sampling technique. 
% * |'wiener-hopf'| selects the method to compute Wiener-Hopf
%   factorizations; can be either |'cr'| (Cyclic Reduction), or |'newton'|
%   (based on the Newton iteration).

%% Syntax
% * |value = cqtoption('option')|
% * |cqtoption('option', value)|

%% Example
%

cqtoption('threshold')