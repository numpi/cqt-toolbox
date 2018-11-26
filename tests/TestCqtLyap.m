function TestCqtLyap
%TESTCQTLYAP Check the implementation of cqtlyap

n=randn(1,6); p=randn(1,4); n(1)=(1 + sum(abs([n,p]))); p(1)=n(1);
T = cqt(n, p, randn(10,2), rand(10,2));
T = T * T';

n=randn(1,6); p=randn(1,4); n(1)=(1 + sum(abs([n,p]))); p(1)=n(1);
S = cqt(n, p, randn(10,2), rand(10,2));
S = S * S';

S = S / norm(S);
T = T / norm(T);

C = cqt([ 1 rand(1,4) ], [1 rand(1,5) ], rand(9, 7));

X = cqtlyap(T, S, C);

CheckTestResult(norm(T*X + X*S + C), '<', ...
    1e4 * cqtoption('threshold') * norm(X), ...
    'CQT Sylvester solver');