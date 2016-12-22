function T2 = conj(T)
%CONJ Conjugate the matrix T. 

if T.sz(1) == inf
    T2 = cqt(conj(T.n), conj(T.p), conj(T.U), conj(T.V));
else
    T2 = cqt(conj(T.n), conj(T.p), conj(T.U), conj(T.V), ...
        conj(T.W), conj(T.Z), T.sz(2), T.sz(1));
end

end

