function TestCqtExpLog
%TESTCQTEXPLOG Check the implementation of expm and logm

n=randn(1,6); p=randn(1,4); n(1)=sum(abs([n,p])); p(1)=n(1);
T = cqt(n, p, randn(10,2), rand(10,2));

% We check the multiplication on a 'large enough' matrix
S = expm(T);

XS = S(1:300, 1:300);
XT = expm(T(1:300,1:300));

res = norm(XS(1:10,1:10) - XT(1:10,1:10));
res2 = norm(T * S + cqt(-1,-1,0,0));
res3 = norm(S * T + cqt(-1,-1,0,0));

CheckTestResult(res, '<', 1e4 * cqtoption('threshold') * norm(S), ...
    'CQT exponential');

[T1] = GenerateFiniteExample(10, 2, 4, 100, 100);
T1 = T1 + cqt(5,5,0,0,0,0,100,100);
T2 = expm(T1);
fT = expm(full(T1));
res = norm(full(T2) - fT) / norm(full(fT));

CheckTestResult(res, '<', 1e4 * cqtoption('threshold'), ...
    'CQT exponential (finite case)');

% We check the multiplication on a 'large enough' matrix
S = logm(T);

XS = S(1:300, 1:300);
XT = logm(T(1:300,1:300));

res = norm(XS(1:10,1:10) - XT(1:10,1:10));

CheckTestResult(res, '<', 1e4 * cqtoption('threshold') * norm(S), ...
    'CQT logarithm');

[T1] = GenerateFiniteExample(10, 2, 4, 100, 100);
T1 = T1 + cqt(5,5,0,0,0,0,100,100);
T2 = logm(T1);
fT = logm(full(T1));
res = norm(full(T2) - fT) / norm(full(fT));

CheckTestResult(res, '<', 1e4 * cqtoption('threshold'), ...
    'CQT logarithm (finite case)');


end
