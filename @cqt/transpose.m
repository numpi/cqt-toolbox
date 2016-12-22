function T2 = transpose(T)
%TRANSPOSE Tranpose the CQT matrix T. 

if T.sz(1) == inf
    T2 = cqt(T.p, T.n,  T.V, T.U);
else
    T2 = cqt(T.p, T.n,  T.V, T.U, T.Z, T.W, T.sz(2), T.sz(1));
end

end

