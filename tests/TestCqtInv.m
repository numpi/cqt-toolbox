function TestCqtMtimes
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
fprintf('TestCqtInv1: Residue on CQT inverse: %e\n', res);
fprintf('TestCqtInv2: Residue on CQT inverse: %e\n', res2);
fprintf('TestCqtInv3: Residue on CQT inverse: %e\n', res3);
res=max([res,res2,res3]);
assert(res < 1000 * eps);

end
