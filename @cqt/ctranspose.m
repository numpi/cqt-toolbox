function T2 = ctranspose(T)
%CTRANSPOSE Conjugate transpose of the matrix T. 

if T.sz(1) == inf
    T2 = cqt(conj(T.p), conj(T.n), conj(T.V), conj(T.U));
else
    T2 = cqt(conj(T.p), conj(T.n), conj(T.V), conj(T.U), ...
        conj(T.Z(end:-1:1,end:-1:1)), ...
        conj(T.W(end:-1:1,end:-1:1)), T.sz(2), T.sz(1));
end


end

