function RunAllTests
%RUNALLTESTS Run all the unit tests. 

addpath ../

% All the tests are repeated with all the possible configurations. At the
% moment we support CR and FFT based inversions. 
inversion = { 'fft', 'cr' };

for i = 1 : length(inversion)
    fprintf('\nCONFIGURATION: Setting inversion = %s\n\n', inversion{i});
    cqtoption('inversion', inversion{i});

    TestCqtGeneric;
    TestCqtTranspose;
    TestCqtPlus;
    TestCqtMtimes;
    TestCqtInv;
    TestCqtUminus;
    TestCqtMpower;
    TestCqtMldivide;
    TestCqtMrdivide;
    TestCqtUL;
    TestCqtEmpty;

end

rmpath ../

end

