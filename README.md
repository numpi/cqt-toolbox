# Quasi Toeplitz matrix toolbox for MATLAB

[![Build status](https://api.travis-ci.org/numpi/cqt-toolbox.svg?branch=master)](https://travis-ci.org/numpi/cqt-toolbox)

The CQT toolbox provides a new MATLAB type that represents semi-infinite
matrices obtained as the sum of a Toeplitz matrix and a finite-dimensional
correction in the top-left corner. 

The arithmetic, plus many additional functions, are implemented (or 
will be included soon) in the toolbox. 

To install the latest version in your MATLAB environment you
can run the following commands: 

```Matlab
unzip('https://github.com/numpi/cqt-toolbox/archive/master.zip'); ... 
movefile('cqt-toolbox-master', 'cqt-toolbox'); ... 
addpath(fullfile('cqt-toolbox')); savepath;
```
   
