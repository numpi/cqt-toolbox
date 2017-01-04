function TestCqtUL
%TESTCQTUL Test the UL factorization

T = GenerateExample(6, 2, 3);

[U,L,E] = ul(T);
CheckTestResult(norm(U*L+E-T), '<', 1e3 * eps * norm(T), ...
    'Accuracy of the UL factorization');

T = GenerateFiniteExample(6, 2, 3, 45, 45);
[U,L,E] = ul(T);
CheckTestResult(norm(U*L+E-T), '<', 1e3 * eps * norm(T, 1), ...
    'Accuracy of the UL factorization (finite case)');

end

