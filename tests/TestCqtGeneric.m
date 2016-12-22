function TestCqtGeneric
%TESTCQTGENERIC Generic tests for the cqt class. 

[T1, ~, U1, V1, sn1, sp1] = GenerateExample(6, 3, 4);

[ssn1, ssp1] = symbol(T1);
res = norm(sn1 - ssn1) + norm(sp1 - ssp1);
fprintf('TestCqtGeneric: Residue on the retrieve symbols: %e\n', res);
assert(res < eps);

[U, V] = correction(T1);
res = norm(U - U1) + norm(V - V1);
fprintf('TestCqtGeneric: Residue on the retrieved factored correction: %e\n', res);
assert(res < eps);

E = correction(T1);
res = norm(E - U1 * V1.');
fprintf('TestCqtGeneric: Residue on the retrieved dense correction: %e\n', res);
assert(res < eps);


end

