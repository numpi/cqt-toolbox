function TestCqtTranspose
%TESTCQTTRANSPOSE Test the implementation of the transposition. 

T = GenerateExample(6, 2, 3);

FT = T(1:100, 1:100);
FTT = T.';
FTT = FTT(1:100, 1:100);

CheckTestResult(norm(FT.' - FTT), '<', 4.0 * eps, ...
    'Transposition operator');

FT = T(1:100, 1:100);
FTT = T';
FTT = FTT(1:100, 1:100);
CheckTestResult(norm(FT' - FTT), '<', 4.0 * eps, ...
    'Conjugate transpose operator');

FT = T(1:100, 1:100);
FTT = conj(T);
FTT = FTT(1:100, 1:100);
CheckTestResult(norm(conj(FT) - FTT), '<', 4.0 * eps, ...
    'Conjugation');

% Test the finite case
T = GenerateFiniteExample(6, 2, 3);

FT = full(T);
FTT = T.';
FTT = full(FTT);
CheckTestResult(norm(FT.' - FTT), '<', 4.0 * eps, ...
    'Transposition operator (finite case)');

FT = full(T);
FTT = T';
FTT = full(FTT);
CheckTestResult(norm(FT' - FTT), '<', 4.0 * eps, ...
    'Conjugate transpose operator (finite case)');

FT = full(T);
FTT = conj(T);
FTT = full(FTT);
CheckTestResult(norm(conj(FT) - FTT), '<', 4.0 * eps, ...
    'Conjugation (finite case)');


end

