function T2 = conj(T)
%CONJ Conjugate the matrix T. 

if T.sz(1) == inf
    T2 = cqt(conj(T.n), conj(T.p), conj(T.U), conj(T.V));
else
    T2 = cqt(conj(T.n), conj(T.p), conj(T.U), conj(T.V), ...
        conj(T.W(end:-1:1,end:-1:1)), ...
        conj(T.Z(end:-1:1,end:-1:1)), T.sz(1), T.sz(2));
end

end

