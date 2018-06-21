function res=qt_norm(am, ap, aU, aV)
% Compute the aqt-norm of a semi-infinite quasi-toeplitz matrix

alfa = ( sqrt(5) + 1 ) / 2;

[~, ru] = qr(aU, 0);
[~ ,rv] = qr(aV, 0);

res = norm(ru * rv');
if ~isempty(am)
    res = res + alfa * sum(abs([ am, ap(2:end) ]));
end

%temp = abs(aU * aV.');
%res = res + sum(sum(temp,1));
