function T2 = ctranspose(T)
%CTRANSPOSE Conjugate transpose of the matrix T. 

if T.sz(1) == inf
    T2 = cqt(conj(T.p), conj(T.n), conj(T.V), conj(T.U));
else
    T2 = cqt(conj(T.p), conj(T.n), conj(T.V), conj(T.U), ...
        conj(T.Z), conj(T.W), T.sz(2), T.sz(1));
end


end

