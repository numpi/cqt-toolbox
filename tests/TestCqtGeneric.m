function TestCqtGeneric
%TESTCQTGENERIC Generic tests for the cqt class. 

[T1, ~, U1, V1, sn1, sp1] = GenerateExample(6, 3, 4);

[ssn1, ssp1] = symbol(T1);
CheckTestResult(norm(sn1 - ssn1) + norm(sp1 - ssp1), '<', eps, ...
    'Accessor method for the symbol');

[U, V] = correction(T1);
CheckTestResult(norm(U - U1) + norm(V - V1), '<', eps, ...
    'Accessor method for the finite correction');

E = correction(T1);
CheckTestResult(norm(E - U1 * V1.'), '<', eps, ...
    'Accessor method for the finite and dense correction');

% Similar tests for the finite case

[T1, ~, U1, V1, W1, Z1, sn1, sp1] = GenerateFiniteExample(6, 3, 4);

[ssn1, ssp1] = symbol(T1);
res = norm(sn1 - ssn1) + norm(sp1 - ssp1);
fprintf('TestCqtGeneric: Residue on the retrieved symbols (finite case): %e\n', res);
assert(res < eps);

[U, V, W, Z] = correction(T1);
res = norm(U - U1) + norm(V - V1) + norm(W - W1) + norm(Z - Z1);
fprintf('TestCqtGeneric: Residue on the retrieved factored correction (finite case): %e\n', res);
assert(res < 4.0 * eps);

[E, F] = correction(T1);
res = norm(E - U1 * V1.') + norm(F - W1 * Z1.');
fprintf('TestCqtGeneric: Residue on the retrieved dense correction: %e\n', res);
assert(res < 4.0 * eps);


end

