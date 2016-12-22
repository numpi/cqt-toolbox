function TestCqtMtimes
%TESTCQTMTIMES Check the implementation of cqt/mtimes

[T1, T2] = GenerateExample(6, 2, 3);

% We check the multiplication on a 'large enough' matrix
S = T1 * T2;

XS = S(1:300, 1:300);
XT1 = T1(1:300,1:300);
XT2 = T2(1:300,1:300);
XT12 = XT1 * XT2;

res = norm(XS(1:10,1:10) - XT12(1:10,1:10));

fprintf('TestCqtMtimes: Residue on CQT multiplication: %e\n', res);
assert(res < 100 * eps);

% TODO: Test the multiplication by scalar.
% Ssr = T1 * 2.0;
% Ssl = 2.0 * T1;




end

