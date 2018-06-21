function TestCqtMldivide
%TESTCQTMLDIVIDE

epsi = cqtoption('threshold');

A = cqt([ 4 -1 ], [4 1], rand(4,2), rand(3,2));
B = cqt(0, 0, rand(4,1), rand(6,1));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * epsi * norm(B) * norm(inv(A)), ...
    'mldivide on CQT and finite matrix');

x = rand;
C = x \ B;
CheckTestResult(norm(C * x - B), '<', 1e2 * epsi * norm(B), ...
    'mldivide with scalar and CQT');

A = cqt([ 4 -1 ], [4 1]);
B = cqt([1 rand(1,2) ], [1 rand(1,2) ]);

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e4 * epsi * norm(inv(A)) * norm(B), ...
    'mldivide on Toeplitz and Toeplitz');

A = cqt([ 4 -1 ], [4 1], randn(4, 2), randn(3, 2));
B = cqt([1 rand(1,2) ], [1 rand(1,3) ], randn(5,3), randn(7, 3));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * epsi * norm(B) * norm(inv(A)), ...
    'mldivide on CQT and CQT');

A = cqt([ 4 -1 ], [4 1], rand(4,2), rand(3,2));
B = cqt(0, 0, rand(4,1), rand(6,1));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * epsi * norm(inv(A)) * norm(B), ...
    'mldivide on CQT and finite matrix');

%
% FINITE CASE
%

m = 100;
n = 125;

A = cqt([ 4 -1 ], [4 1], 0, 0, m, m);
B = cqt([1 rand(1,2) ], [1 rand(1,2) ], 0, 0, m, n);

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * epsi, 'mldivide on Toeplitz and Toeplitz (finite case)');

A = cqt([ 4 -1 ], [4 1], randn(4, 2), randn(3, 2), 0, 0, m, m);
B = cqt([1 rand(1,2) ], [1 rand(1,3) ], randn(5,3), randn(7, 3), 0, 0, m, n);

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * epsi, 'mldivide on CQT and CQT (finite case)');


end

