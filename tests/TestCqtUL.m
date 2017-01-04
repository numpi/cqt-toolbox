function TestCqtUL
%TESTCQTUL Test the UL factorization

T = GenerateExample(6, 2, 3);

[U,L,E] = ul(T);
CheckTestResult(norm(U*L+E-T), '<', eps * 1e3, ...
    'Accuracy of the UL factorization');

T = GenerateFiniteExample(6, 2 ,3, 45, 60);
[U,L,E] = ul(T);
CheckTestResult(norm(U*L+E-T), '<', eps * 1e5, ...
    'Accuracy of the UL factorization (finite case)');

end

