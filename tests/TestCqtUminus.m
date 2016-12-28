function TestCqtUminus
%TESTCQTUMINUS Test the uminus operator. 

T = GenerateExample(6, 2, 3);
FT = - T(1:100,1:100);
T = - T;
res = norm(FT - T(1:100,1:100));

fprintf('TestCqtUminus: Residue on uminus: %e\n', res);
assert(res < eps);

T = GenerateFiniteExample(6, 2, 3);
FT = - full(T);
T = - T;
res = norm(FT - full(T));

fprintf('TestCqtUminus: Residue on uminus (finite case): %e\n', res);
assert(res < eps);


end

