function [ T ] = mtimes(T1, T2)
%MTIMES Multiply two CQT matrices T1 and T2. 
%
%     T = MTIMES(T1, T2) computes the CQT matrix T obtained multiplying two CQT
%     matrices T1 and T2. If T2 is a dense finite matrix it is considered
%     as a semiinfinite matrix with only the leading top-left corner
%     different from zero. 
if isa(T1,'cqt')
	if isa(T2, 'cqt')
    		[cm,cp,cU,cV]=qt_mult(T1.n, T1.p, T1.U, T1.V, T2.n, T2.p, T2.U, T2.V);
    		T = cqt(cm, cp, cU, cV);
	elseif isscalar(T2)
    		T = cqt(T1.n * T2, T1.p * T2, T1.U * T2, T1.V);
	else
		error('Incompatible types multiplication. \nIf you want to multiply a cqt matrix T with a finite matrix of %s A you can use T * cqt(A) ',class(T2));
	end
else
	if isscalar(T1) && isa(T2, 'cqt')
		T = cqt(T2.n * T1, T2.p * T1, T2.U * T1, T2.V);
	else
		error('Incompatible types multiplication. \nIf you want to multiply a finite matrix of %s A with a cqt matrix you can use cqt(A) * T',class(T1));
	end
end

