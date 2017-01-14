function T2 = transpose(T)
%TRANSPOSE Tranpose the CQT matrix T.

T2 = cqt(T.p, T.n, T.V, T.U, T.Z(end:-1:1,end:-1:1), ...
	T.W(end:-1:1,end:-1:1), T.sz(2), T.sz(1));

end

