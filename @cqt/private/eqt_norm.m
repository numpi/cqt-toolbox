function r = eqt_norm(T)
%EQT_NORM Extended QT norm.

[n, p] = symbol(T);
E = correction(T);
e = ones(size(E,1), 1);

E = formatted_sum(E, e * T.c);

if isempty(p)
    rs = 0.0;
elseif length(p) == 1
    rs = norm(n, 1);
else
    rs = norm([ n, p(2:end) ], 1);
end

r = rs + max(norm(E, inf), norm(T.c, inf));

end

