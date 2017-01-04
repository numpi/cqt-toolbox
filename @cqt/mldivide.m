function X = mldivide(A,B)
if isa(A,'cqt') && isa(B,'cqt')
  	[cm, cp, cU, cV] = qt_mldivide(A.n, A.p, A.U, A.V, B.n, B.p, B.U, B.V);
	X = cqt(cm, cp, cU, cV); 
else
	error('A and B must be cqt matrices');
end
