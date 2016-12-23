function disp(T)
%DISP Display a CQT matrix on the screen. 

fprintf('CQT Matrix of size %d x %d\n\n', T.sz(1), T.sz(2));

% Check if the CQT has a non-zero Toeplitz part
if length(T.n) + length(T.p) > 0
    fprintf(' - Toeplitz part (leading %d x %d block): \n', ...
        length(T.n) + 2, length(T.p) + 2);
    disp(toeplitz([ T.n 0 0], [ T.p 0 0 ]));
end

if size(T.U, 1) > 0 && size(T.V, 2) > 0
    fprintf('\n - Finite correction (top-left corner): \n');
    disp(T.U * T.V.');
end

if size(T.W, 1) > 0 && size(T.Z, 2) > 0
    fprintf('\n - Finite correction (bottom-right corner): \n');
    S = T.W * T.Z';
    disp(S(end:-1:1,end:-1:1));
end

end

