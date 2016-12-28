function RunAllTests
%RUNALLTESTS Run all the unit tests. 

addpath ../

TestCqtGeneric;
TestCqtTranspose;
TestCqtPlus;
TestCqtMtimes;
TestCqtInv;
TestCqtUminus;
TestCqtMpower;

rmpath ../

end

