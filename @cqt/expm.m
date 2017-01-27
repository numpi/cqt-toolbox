function eT = expm(T, method)

if max(T.sz) == inf
	nrm_type = 'cqt';
else
	nrm_type = inf;
end

nrm = norm(T, nrm_type);

h = floor(log2(nrm)) + 2;
T = T/2^h;
if ~exist('method','var')
	method = 'pade';
end
if strcmp(method,'taylor')
	maxit = 26;
	eT = cqt(1,1,[],[],T.sz(1),T.sz(2));
	tempT = eT;
	for i=1:maxit
		tempT = tempT * T/i;
		eT = eT + tempT;
	end
elseif strcmp(method,'pade')
	
	c = 1 / 2;
	eTn = cqt(1, 1, [], [], T.sz(1), T.sz(2)) + c*T;
	eTd = cqt(1, 1, [], [], T.sz(1), T.sz(2)) - c*T;
	
	q = 6;
	p = 1;
	X = T;
	for k = 2 : q
		c = c * (q-k+1) / (k*(2*q-k+1));
		X = T * X;
		cX = c*X;
		eTn = eTn + cX;
		if p
			eTd = eTd + cX;
		else
			eTd = eTd - cX;
		end
		p = ~p;
	end
	
	eT =  eTn * inv(eTd);
else
	error('Invalid parameter method in EXPM');
end
eT = eT^(2^h);
