function TestCqtMtimes
%TESTCQTMTIMES Check the implementation of cqt/mtimes

[T1, T2, U1, V1, sn1, sp1, U2, V2, sn2, sp2] = GenerateExample(6, 2, 3);

% We check the multiplication on a 'large enough' matrix
S = T1 * T2;

XS = S(1:300, 1:300);
XT1 = T1(1:300,1:300);
XT2 = T2(1:300,1:300);
XT12 = XT1 * XT2;

res = norm(XS(1:10,1:10) - XT12(1:10,1:10));

fprintf('TestCqtMtimes: Residue on CQT multiplication: %e\n', res);
assert(res < 100 * eps);

Ssr = T1 * 2.0;
Ssl = 2.0 * T1;

[Ssr_n, Ssr_p] = symbol(Ssr);
[Ssl_n, Ssl_p] = symbol(Ssl);

res = norm(Ssr_n - 2.0 * sn1) + norm(Ssr_p - 2.0 * sp1);
fprintf('TestCqtMtimes: Residue on CQT * scalar multiplication: %e\n', ...
    res);
assert(res < eps);

res = norm(Ssl_n - 2.0 * sn1) + norm(Ssl_p - 2.0 * sp1);
fprintf('TestCqtMtimes: Residue on scalar * CQT multiplication: %e\n', ...
    res);
assert(res < eps);

[Ul, Vl] = correction(Ssr);
res = norm(2.0 * U1 * V1' - Ul * Vl');
fprintf('TestCqtMtimes: Residue on CQT * scalar multiplication: %e\n', ...
   res);
assert(res < norm(correction(Ssr) * eps));

[Ul, Vl] = correction(Ssl);
res = norm(2.0 * U1 * V1' - Ul * Vl');
fprintf('TestCqtMtimes: Residue on scalar * CQT multiplication: %e\n', ...
   res);
assert(res < norm(correction(Ssl) * eps));

%% Testing the finite case
[T1, T2, U1, V1, W1, Z1, sn1, sp1, U2, V2, W2, Z2 sn2, sp2] = GenerateFiniteExample(80, 3, 2);
S = T1 * T2;
[v1,v2]=symbol(S);
[E, F] = correction(S);
XT1 = full(T1);
XT2 =full(T2);
XS = full(S);
XT12 = XT1 * XT2;
[v1,v2]=symbol(S);
res = norm(XS - XT12);

%XT12-toep(v1,v2,100,80);
%F
fprintf('TestCqtMtimes: Residue on finite CQT multiplication: %e\n', res);
assert(res < 1000 * eps);
Ssr = T1 * 2.0;
Ssl = 2.0 * T1;

[Ssr_n, Ssr_p] = symbol(Ssr);
[Ssl_n, Ssl_p] = symbol(Ssl);

res = norm(Ssr_n - 2.0 * sn1) + norm(Ssr_p - 2.0 * sp1);
fprintf('TestCqtMtimes: Residue on finite CQT * scalar multiplication: %e\n', ...
    res);
assert(res < eps);
end

