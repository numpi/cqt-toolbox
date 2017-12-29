function [TU, TV] = compress_qr(U, V, nrm)
%[size(U),size(V)]
%COMPRESS_QR Compress the factors of an outer product X1 * X2'.
%
% [Y1, Y2] = COMPRESS_QR(X1, X2) compute a new factorization Y1 * Y2' of
%     the outer product X1 * X2' with (possibly) less columns. It uses a
%     reduced singular value decomposition to compute this optimal
%     representation.
%
% Date: June 2 2016
% Author: Dario A. Bini

if exist('nrm', 'var')
    epsi= nrm * cqtoption('threshold');
else
    epsi = inf;
end

pivoting = 0;
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
if pivoting
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
else
	[q1,r1] = qr(U,0);
	[q2,r2] = qr(V,0);
	r = r1*r2.';        
    
    if epsi == inf
        epsi = norm(r) * cqtoption('threshold');
    end
    
	n1 = size(r1,1);
	n2 = size(r2,1);
end
dosvd = 1;
if dosvd
	[ru,rs,rv] = svd(r);
	
	rs=diag(rs);
	nsv = sum(rs > epsi);
	
	% symmetric scaling
	srs = sqrt(rs(1:nsv));
	TU = q1(:,1:n1)*(ru(:,1:nsv));
	TV = q2(:,1:n2)*(conj(rv(:,1:nsv)));
    
    if exist('nrm', 'var')
        [TU, TV] = svd_clean(TU, TV, rs(1:nsv), nrm);
    end
    TU = TU * diag(srs);
    TV = TV * diag(srs);
	
else        % no svd
	if n2>n1
		TU = q1(:,1:n1);
		TV = q2(:,1:n2)*r';
	else
		TU = q1(:,1:n1)*r;
		TV = q2(:,1:n2);
	end
end

