function disp(T)
%DISP Display a CQT matrix on the screen. 

fprintf('CQT Matrix\n\n');
fprintf(' - Toeplitz part: \n\n');
disp(toeplitz([ T.n 0 0], [ T.p 0 0 ]));
fprintf('\n - Finite part (top-left corner): \n\n');
disp(T.U * T.V.');

end

