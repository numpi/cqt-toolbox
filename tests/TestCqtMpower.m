function TestCqtMpower
%TESTCQTMPOWER Check matrix power of CQT objects. 

T = GenerateExample(6, 2, 3);

S = T(1:1000, 1:1000);
T2 = T^2;
S2 = S^2;

CheckTestResult(norm(T2(1:100,1:100) - S2(1:100,1:100)), '<', ...
    1e2 * eps * norm(S2(1:100,1:100)), ...
    'Computation of T^2 for infinite CQT matrices');

T3 = T^3;
S3 = S^3;

CheckTestResult(norm(T3(1:100,1:100) - S3(1:100,1:100)), '<', ...
    1e2 * eps * norm(S3(1:100,1:100)), ...
    'Computation of T^3 for infinite CQT matrices');

T13 = T^5;
S13 = S^5;

CheckTestResult(norm(T13(1:100,1:100) - S13(1:100,1:100)), '<', ...
    1e2 * eps * norm(S13(1:100,1:100)), ...
    'Computation of T^5 for infinite CQT matrices');

T = T + cqt(5, 5);

S = T(1:1000, 1:1000);
T2 = T^(-2);
S2 = S^(-2);

CheckTestResult(norm(T2(1:100,1:100) - S2(1:100,1:100)), '<', ...
    1e2 * eps * norm(S2(1:100,1:100)) + norm(T) * eps, ...
    'Computation of T^(-2) for infinite CQT matrices');

T3 = T^(-3);
S3 = S^(-3);

CheckTestResult(norm(T3(1:100,1:100) - S3(1:100,1:100)), '<', ...
    1e2 * eps * norm(S3(1:100,1:100)) + norm(T) * eps, ...
    'Computation of T^(-3) for infinite CQT matrices');

T13 = T^(-5);
S13 = S^(-5);

CheckTestResult(norm(T13(1:100,1:100) - S13(1:100,1:100)), '<', ...
    1e2 * eps * norm(S13(1:100,1:100)) + norm(T) * eps, ...
    'Computation of T^(-5) for infinite CQT matrices');

%
% Finite case
%

T = GenerateFiniteExample(6, 2, 3, 100, 100);
T = T + cqt(5, 5, 0, 0, 0, 0, 100, 100);
S = full(T);

for p = [ 1, 2, 3, 5, -1, -2, -3, -5 ]
    TT = T^p;
    SS = S^p;

    CheckTestResult(norm(full(TT) - SS), '<', 1e3 * eps * norm(SS), ...
       sprintf('Computation of T^(%d) for finite CQT matrices', p));
end

end

