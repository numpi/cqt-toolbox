function [Am1, A0, A1, hA0] = jackson(N)

switch N 
	case 1
		l1=1;l2=0;mu1=1.5;mu2=2;p=1;q=0;
	case 2
		l1=1;l2=0;mu1=2;mu2=1.5;p=1;q=0;
	case 3
		l1=0;l2=1;mu1=1.5;mu2=2;p=0;q=1;
	case 4
		l1=0;l2=1;mu1=2;mu2=1.5;p=0;q=1;
	case 5
		l1=1;l2=1;mu1=2;mu2=2;p=0.1;q=0.8;
	case 6
		l1=1;l2=1;mu1=2;mu2=2;p=0.8;q=0.1;
	case 7
		l1=1;l2=1;mu1=2;mu2=2;p=0.4;q=0.4;
	case 8
		l1=1;l2=1;mu1=10;mu2=10;p=0.5;q=0.5;
	case 9
		l1=1;l2=5;mu1=10;mu2=15;p=0.4;q=0.9;
	case 10
		l1=5;l2=1;mu1=15;mu2=10;p=0.9;q=0.4;
end

A0 = cqt([ l1 + l2 + mu1 + mu2, (p - 1) * mu1 ], ...
	[ l1 + l2 + mu1 + mu2, -l1 ], ...
	-mu1);

A1 = cqt([ -l2 , -p*mu1 ], -l2);

Am1 = cqt((q-1)*mu2, [(q-1)*mu2, -q*mu2]);

hA0 = cqt([ l1 + l2 + mu1, (p-1)*mu1 ], ...
	[ l1 + l2 + mu1, -l1 ], -mu1);

% Perform normalization of the rows 
f  = 1 / (l1 + l2 + mu1 + mu2);
A0 = f * A0;
A1 = f * A1;
Am1 = f * Am1;
hA0 = f * hA0;
