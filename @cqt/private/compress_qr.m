function [TU, TV] = compress_qr(U, V)

%COMPRESS_QR Compress the factors of an outer product X1 * X2'.
%
% [Y1, Y2] = COMPRESS_QR(X1, X2) compute a new factorization Y1 * Y2' of
%     the outer product X1 * X2' with (possibly) less columns. It uses a
%     reduced singular value decomposition to compute this optimal
%     representation.
%
% Date: June 2 2016
% Author: Dario A. Bini
epsi= eps;

if size(U,1) == 1
	TU = 1;  TV = V*U.';
	return
end

if isempty(U) || isempty(V)
	TU = [];
	TV = [];
	return;
end

if max(max(abs(U)))==0 || max(max(abs(V)))==0
	TU = [];  TV = [];
	return
end
[q1,r1,p1] = qr(U,0);
[q2,r2,p2] = qr(V,0);

% invert the permutations
for i=1:length(p1)
	ip1(p1(i)) = i;
end
for i=1:length(p2)
	ip2(p2(i)) = i;
end

% compute the number of meaningful elements
% gestire il caso di r1 vuoto
n1 = sum(abs(diag(r1))>epsi*abs(r1(1,1)));
n2 = sum(abs(diag(r2))>epsi*abs(r2(1,1)));
r = r1(1:n1,ip1)*r2(1:n2,ip2).';
if isempty(r)
	TU=[];
	TV=[];
	return;
end
dosvd = 1;
if dosvd
	[ru,rs,rv] = svd(r);
	
	rs=diag(rs);
	nsv = sum(rs>epsi*rs(1));
	
	% symmetric scaling
	srs = sqrt(rs(1:nsv));
	TU = q1(:,1:n1)*(ru(:,1:nsv)*diag(srs));
	TV = q2(:,1:n2)*(conj(rv(:,1:nsv))*diag(srs));
	
else        % no svd
	if n2>n1
		TU = q1(:,1:n1);
		TV = q2(:,1:n2)*r';
	else
		TU = q1(:,1:n1)*r;
		TV = q2(:,1:n2);
	end
end
e1 = abs(TU)*ones(size(TU,2),1);
e2 = abs(TV)*ones(size(TV,2),1);
e1 = e1>eps*max(e1);
e2 = e2>eps*max(e2);
n1 = size(TU,1);
for i=length(e1):-1:1
	if e1(i)
		break
	else
		n1 = n1-1;
	end
end
n2 = size(TV,1);
for i=length(e2):-1:1
	if e2(i)
		break
	else
		n2 = n2-1;
	end
end
TU = TU(1:n1,:); TV = TV(1:n2,:);


