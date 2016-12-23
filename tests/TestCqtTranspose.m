function TestCqtTranspose
%TESTCQTTRANSPOSE Test the implementation of the transposition. 

T = GenerateExample(6, 2, 3);

FT = T(1:100, 1:100);
FTT = T.';
FTT = FTT(1:100, 1:100);
res = norm(FT.' - FTT);

fprintf('TestCqtTranspose: Residue on transposition: %e\n', res);
assert(res < 4.0 * eps);

FT = T(1:100, 1:100);
FTT = T';
FTT = FTT(1:100, 1:100);
res = norm(FT' - FTT);

fprintf('TestCqtTranspose: Residue on (conjugate) transposition: %e\n', res);
assert(res < 4.0 * eps);

FT = T(1:100, 1:100);
FTT = conj(T);
FTT = FTT(1:100, 1:100);
res = norm(conj(FT) - FTT);

fprintf('TestCqtTranspose: Residue on conjugation: %e\n', res);
assert(res < 4.0 * eps);

% Test the finite case
T = GenerateFiniteExample(6, 2, 3);

FT = full(T);
FTT = T.';
FTT = full(FTT);
res = norm(FT.' - FTT);

fprintf('TestCqtTranspose: Residue on transposition (finite case): %e\n', res);
assert(res < 4.0 * eps);

FT = full(T);
FTT = T';
FTT = full(FTT);
res = norm(FT' - FTT);

fprintf('TestCqtTranspose: Residue on (conjugate) transposition (finite case): %e\n', res);
assert(res < 4.0 * eps);

FT = full(T);
FTT = conj(T);
FTT = full(FTT);
res = norm(conj(FT) - FTT);

fprintf('TestCqtTranspose: Residue on conjugation (finite case): %e\n', res);
assert(res < 4.0 * eps);


end

