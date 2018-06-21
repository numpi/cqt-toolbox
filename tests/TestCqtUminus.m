function TestCqtUminus
%TESTCQTUMINUS Test the uminus operator.

T = GenerateExample(6, 2, 3);
FT = - T(1:100,1:100);
T = - T;
res = norm(FT - T(1:100,1:100));

CheckTestResult(res, '<', eps, ...
    'Uminus implementation');

T = GenerateFiniteExample(6, 2, 3);
FT = - full(T);
T = - T;
res = norm(FT - full(T));

CheckTestResult(res, '<', eps, ...
    'Uminus implementation (finite case)');

end

