function RunAllTests
%RUNALLTESTS Run all the unit tests. 

addpath ../

TestCqtGeneric;
TestCqtPlus;
TestCqtMtimes;
TestCqtInv;

rmpath ../

end

