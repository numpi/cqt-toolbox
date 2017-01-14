function TestCqtEmpty
%TESTCQTEMPTY Test the arithmetic of empty matrices. 

A = cqt([], [], [], [], inf, 0);
B = cqt([], [], [], [], 0, inf);

C = A * B;

CheckTestResult(size(A, 1), '==', inf, ...
    'Size on empty matrices (inf x 0), left dimension');
CheckTestResult(size(A, 2), '==', 0,   ...
    'Size on empty matrices (inf x 0), right dimension');
CheckTestResult(size(B, 1), '==', 0,   ...
    'Size on empty matrices (0 x inf), left dimension');
CheckTestResult(size(B, 2), '==', inf, ...
    'Size on empty matrices (0 x inf), right dimension');
CheckTestResult(size(C, 1), '==', inf, ...
    'Size on empty matrices (inf x 0) * (0 x inf), left dimension');
CheckTestResult(size(C, 2), '==', inf, ...
    'Size on empty matrices (inf x 0) * (0 x inf), right dimension');
CheckTestResult(norm(C), '==', 0, ...
    'Norm of outer product of empty matrices == 0')
CheckTestResult(isempty(full(B * A)), '==', true, ...
    'Inner product of empty matrices == []')


end

