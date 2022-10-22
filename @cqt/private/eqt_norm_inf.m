function r = eqt_norm_inf(T)
%EQT_NORM_INF Infinity norm of the extended QT matrix T. 

[am, ap] = symbol(T);

if ~isempty(am)
    am = am(2:end);
end

r = norm([ am, ap ], 1) + norm(T.c, 1);

% Find the number of rows where correction and the Toeplitz part overlap
kr = max(size(T.V, 1), length(T.c)) + length(am);

% Take the maximum with the norm of the rows in the first rows -- only goes
% on the support of the correction. 
bs = 128;
for i = 1 : bs : min( size(T.U, 1), kr )
    j = max([ length(ap) + i, length(T.c), size(T.V, 1) ]);    
    RR = slicemat(T, { i : i+bs-1 , 1 : j });
    r = max(r, norm(RR, inf));
end

% Consider the remaining part of the correction, if any
r = max(r, ...
    norm(formatted_sum(...
        T.U(kr+1:end, :) * T.V.', ...
        ones(size(T.U, 1) - kr, 1) * T.c), inf) ...
    + norm([ am, ap ], 1));


end

