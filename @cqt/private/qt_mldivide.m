function C = qt_mldivide(A, B)
%QT_MLDIVIDE Compute inv(A) * B.

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
    C = inv(A) * B;
    return;
end

[U, L, E] = ul(A);

% Factor E as the outer product E1 * E2.', with E1 and E2 CQT matrices.
E1 = cqt([], [], E.U, [], inf, size(E.U, 2));
E2 = cqt([], [], E.V, [], inf, size(E.V, 2));

% Compute A \ B by Sherman-Morrison

Linv = inv(L);
Uinv = inv(U);

explicit_inverse = false;

if explicit_inverse
    LUinv = Linv * Uinv;
    LUB  = LUinv * B;
    LUE1 = LUinv * E1;
else
    LUB = Linv * (Uinv * B);    
    
    % LUE1 = Linv * (Uinv * E1);
    E1U = E.U; nrmE = norm(E1U);
    E1U = toepmult_fft(Uinv.n, Uinv.p, ...
          size(E1U, 1) + length(Uinv.p) - 1, size(E1U, 1), E1U);
    [E1U, ~] = svd_clean(E1U, eye(size(E1U, 2)), eye(size(E1U, 2)), nrmE * norm(Uinv.p, 1));
    E1U = toepmult_fft(Linv.n, Linv.p, ...
        size(E1U, 1) + length(Linv.n) - 1, size(E1U, 1), E1U);
    LUE1 = cqt([], [], E1U, [], inf, size(E1U, 2));  
end

% S = eye(size(E.U, 2)) + full(E2.' * LUE1);
% E1U = LUE1.U; 
if size(E1U, 1) < size(E.V, 1)
    E1U(size(E.V, 1), 1) = 0;
else
    E1U = E1U(1:size(E.V, 1), :);
end
S = eye(size(E.V, 2)) + ( (E.V).' * E1U ) * LUE1.V.';

% LUE1SE2 = - LUE1 * (S \ E2.');
SV = LUE1.V.' * (S \ E.V.');
SU = -LUE1.U;
LUE1SE2 = cqt([], [], SU, SV.');

LUE1SE2.p = 1;
LUE1SE2.n = 1;
C = LUE1SE2 * LUB;
