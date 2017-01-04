function res=qt_norm(am, ap, aU, aV)
% Compute the aqt-norm of a semi-infinite quasi-toeplitz matrix
res = 0;
for j=1:length(am)
	res = res + j* abs(am(j));
end
for j=2:length(ap)
	res = res + j * abs(ap(j));
end

temp = abs(aU * aV.');
res = res + sum(sum(temp,1));
