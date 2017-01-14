function [l, u] = spectral_cr(vm, vp)
%function [l,u]=spectral_cr(vm,vp)
% Computes the spectral factorization of the Laurent
% polynomial with coefficients vm, vp, by means of CR

% Handle the trivial case of a diagonal Toeplitz matrix
if (length(vm) == 1) && (length(vp) == 1)
	l = 1;
	u = vp;
	return;
end

maxiter = 40;  epsi = 1.e-20;
if size(vm,1)==1
	vm = vm.';
end
if size(vp,1)==1
	vp = vp.';
end
nm = length(vm); np = length(vp);
if nm>np
	am = vm(1:nm-1); ap = [vp;zeros(nm-np-1,1)];
	A = toeplitz(am,ap);
	B = toeplitz([vm(nm);zeros(nm-2,1)],vm(nm:-1:2));
	C = toeplitz([zeros(nm-np,1);vp(np:-1:2)],zeros(nm-1,1));
elseif nm<np
	am = [vm;zeros(np-nm-1,1)]; ap = vp(1:np-1);
	A = toeplitz(am,ap);
	B = toeplitz(zeros(np-1,1),[zeros(np-nm,1);vm(nm:-1:2)]);
	C = toeplitz(vp(np:-1:2),[vp(np);zeros(np-2,1)]);
else
	am = vm(1:nm-1); ap = [vp(1:np-1)];
	A = toeplitz(am,ap);
	B = toeplitz([vm(nm);zeros(nm-2,1)],vm(nm:-1:2));
	C = toeplitz(vp(np:-1:2),[vp(np);zeros(np-2,1)]);
end

%CR
b = B; c = C;
n = size(A,1);
at = A;  ah = A;
for k=1:maxiter
	ABC = A\[B,C];
	cab = C*ABC(:,1:n);
	bac = B*ABC(:,n+1:2*n);
	A = A-cab-bac;
	B = -B*ABC(:,1:n);
	C = -C*ABC(:,n+1:2*n);
	at = at-cab;
	ah = ah-bac;
	err = sqrt(norm(B,'inf')*norm(C,'inf'));
	% err
	if err<epsi
		break
	end
end

en = zeros(n,1); e1 = en; en(n) = 1; e1(1) = 1;
ll = b*(ah\en);
l(1) = 1; l(2:n+1) = ll; l = l.';
u = c*(at\e1);
u(n+1) = 1;
l = l(1:nm);
u = u(end:-1:end-np+1);

m = min(length(l), length(u));
y = vp(1) / dot(l(1:m), u(1:m));
x = sqrt(abs(y));
u = u*x*sign(y); l = l*x;

% inv(A)
end

