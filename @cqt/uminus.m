function T = uminus(T1)
%UMINUS change the sign of a CQT-matrix. 
%
%     T = UMINUS(T1) changes the sign of the CQT-matrix T1 and produces a new CQT matrix T. 

    T = cqt(-T1.n, -T1.p, -T1.U, T1.V);

end
