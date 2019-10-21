function r = eqt_norm_inf(T)
%EQT_NORM_INF Infinity norm of the extended QT matrix T. 

[am, ap] = symbol(T);

if ~isempty(am)
    am = am(2:end);
end

r = norm([ am, ap ], 1) + norm(T.c, 1);

% Take the maximum with the norm of the rows in the first rows -- only goes
% on the support of the correction. 
bs = 128;
for i = 1 : bs : size(T.U, 1)
    j = max([ length(ap) + i, length(T.c), size(T.V, 1) ]);    
    RR = slicemat(T, { i : i+bs-1 , 1 : j });
    r = max(r, norm(RR, inf));
end


end

