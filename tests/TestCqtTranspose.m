function TestCqtTranspose
%TESTCQTTRANSPOSE Test the implementation of the transposition. 

T = GenerateExample(6, 2, 3);

FT = T(1:100, 1:100);
FTT = T.';
FTT = FTT(1:100, 1:100);
res = norm(FT.' - FTT);

fprintf('Residue on transposition: %e\n', res);
assert(res < eps);

FT = T(1:100, 1:100);
FTT = T';
FTT = FTT(1:100, 1:100);
res = norm(FT' - FTT);

fprintf('Residue on (conjugate) transposition: %e\n', res);
assert(res < eps);

FT = T(1:100, 1:100);
FTT = conj(T);
FTT = FTT(1:100, 1:100);
res = norm(conj(FT) - FTT);

fprintf('Residue on conjugation: %e\n', res);
assert(res < eps);


% Test the finite case


end

