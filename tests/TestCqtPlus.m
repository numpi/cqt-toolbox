function TestCqtPlus
%TESTCQTPLUS Check the implementation of cqt/plus
%
% Author: Leonardo Robol <leonardo.robol@cs.kuleuven.be>

epsi = cqtoption('threshold');

% Size of the top-left correction;
[T1, T2] = GenerateExample(6, 2, 3);

S = T1 + T2;

res = norm(full(S(1:20,1:20)) - full(T1(1:20,1:20)) - full(T2(1:20,1:20)));

CheckTestResult(res, '<', 1e3 * epsi * norm(S), ...
    'Sum of CQT matrices');

[T1, T2] = GenerateFiniteExample(6, 2, 3);
S = T1 + T2;
res = norm(full(S) - full(T1) - full(T2));

CheckTestResult(res, '<', 1e3 * epsi * norm(S), ...
    'Sum of CQT matrices (finite case)');

end

