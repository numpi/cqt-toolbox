function [k, v] = epslength(a)
%EPSLENGTH Number of non-negligible entries in a.
%
% determine the numerical length of the vector and remove the
% trailing entries if sum(a(k+1:end))<eps*norm(a,'inf')
%  global epsi
%  global relative

epsi = eps;
relative = true;
if relative
	mx=norm(a,'inf');
	epsx=mx*epsi;
else
	epsx=epsi;
end
b=a(end:-1:1);
n=length(b);
s=0;
for i=1:n
	s=s+abs(b(i));
	if s>epsx
		break
	end
end
k=n-i+1;
v=a(1:k);
