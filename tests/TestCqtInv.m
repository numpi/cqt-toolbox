function TestCqtInv
%TESTCQTINV Check the implementation of cqt/inv

n=randn(1,6); p=randn(1,4); n(1)=sum(abs([n,p])); p(1)=n(1);
T = cqt(n, p, randn(10,2), rand(10,2));

% We check the multiplication on a 'large enough' matrix
S = inv(T);

XS = S(1:300, 1:300);
XT = inv(T(1:300,1:300));

res = norm(XS(1:10,1:10) - XT(1:10,1:10));
res2 = norm(T * S + cqt(-1,-1,0,0));
res3 = norm(S * T + cqt(-1,-1,0,0));

CheckTestResult(res, '<', 1e4 * cqtoption('threshold') * norm(S), ...
	'CQT inversion');
CheckTestResult(res2, '<', 1e4 * cqtoption('threshold') * norm(S), ...
	'CQT inversion');
CheckTestResult(res3, '<', 1e4 * cqtoption('threshold') * norm(S), ...
	'CQT inversion');

[T1] = GenerateFiniteExample(10, 2, 4, 100, 100);
T1 = T1 + cqt(5,5,0,0,0,0,100,100);
T2 = inv(T1);
fT = inv(full(T1));
res = norm(full(T2) - fT);
% fprintf('TestCqtInv4: Residue on a finite CQT inverse: %e\n', res3);
% assert(res < 1e4 * eps);

CheckTestResult(res, '<', 1e4 * cqtoption('threshold'), ...
	'CQT inversion (finite case)');


end
