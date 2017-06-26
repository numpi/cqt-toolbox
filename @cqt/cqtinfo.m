function cqtinfo(C)
%CQTINFO Get information on a CQT object. 

fprintf('\n CQT matrix of size %d x %d\n', size(C, 1), size(C, 2));
fprintf('  - Rank of the correction: %d\n', cqtrank(C));
fprintf('  - Length of positive / negative symbol: %d / %d\n', ...
    length(C.p), length(C.n));

fprintf('\n')

end

