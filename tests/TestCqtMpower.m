function TestCqtMpower
%TESTCQTMPOWER Check matrix power of CQT objects. 

T = GenerateExample(6, 2, 3);

S = T(1:1000, 1:1000);
T2 = T^2;
S2 = S^2;

CheckTestResult(norm(T2(1:100,1:100) - S2(1:100,1:100)), '<', 1e4 * eps, ...
    'Computation of T^2 for infinite CQT matrices');

T3 = T^3;
S3 = S^3;

CheckTestResult(norm(T3(1:100,1:100) - S3(1:100,1:100)), '<', 1e4 * eps, ...
    'Computation of T^3 for infinite CQT matrices');

T13 = T^13;
S13 = S^13;

CheckTestResult(norm(T13(1:100,1:100) - S13(1:100,1:100)), '<', 1e-5, ...
    'Computation of T^13 for infinite CQT matrices');

% The finite multiplication has not been implemented, yet. 
return;

%
% Finite case
%

T = GenerateFiniteExample(6, 2, 3, 100, 100);

S = full(T);
T2 = T^2;
S2 = S^2;

CheckTestResult(norm(full(T2) - S2), '<', 1e4 * eps, ...
    'Computation of T^2 for infinite CQT matrices');

T3 = T^3;
S3 = S^3;

CheckTestResult(norm(full(T3) - S3), '<', 1e4 * eps, ...
    'Computation of T^3 for infinite CQT matrices');

T13 = T^13;
S13 = S^13;

CheckTestResult(norm(full(T13) - S13), '<', 1e-6, ...
    'Computation of T^13 for infinite CQT matrices');

end

