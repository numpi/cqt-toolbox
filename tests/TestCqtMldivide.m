function TestCqtMldivide
%TESTCQTMLDIVIDE

A = cqt([ 4 -1 ], [4 1], rand(4,2), rand(3,2));
B = cqt(0, 0, rand(4,1), rand(6,1));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * eps * norm(B) * norm(inv(A)), ...
	'mldivide on CQT and finite matrix');

A = cqt([ 4 -1 ], [4 1]);
B = cqt([1 rand(1,2) ], [1 rand(1,2) ]);

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * eps * norm(inv(A)) * norm(B), ...
	'mldivide on Toeplitz and Toeplitz');

A = cqt([ 4 -1 ], [4 1], randn(4, 2), randn(3, 2));
B = cqt([1 rand(1,2) ], [1 rand(1,3) ], randn(5,3), randn(7, 3));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * eps * norm(B) * norm(inv(A)), ...
	'mldivide on CQT and CQT');

A = cqt([ 4 -1 ], [4 1], rand(4,2), rand(3,2));
B = cqt(0, 0, rand(4,1), rand(6,1));

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * eps * norm(inv(A)) * norm(B), ...
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

CheckTestResult(norm(C - D), '<', 1e3 * eps, 'mldivide on Toeplitz and Toeplitz (finite case)');

A = cqt([ 4 -1 ], [4 1], randn(4, 2), randn(3, 2), 0, 0, m, m);
B = cqt([1 rand(1,2) ], [1 rand(1,3) ], randn(5,3), randn(7, 3), 0, 0, m, n);

C = inv(A) * B;
D = A \ B;

CheckTestResult(norm(C - D), '<', 1e3 * eps, 'mldivide on CQT and CQT (finite case)');


end

