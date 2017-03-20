% This script applies Cyclic Reduction to solve the matrix equation
%                 B*G^2 + A*G + C=0
% for some test problems where the blocks A, B, C, are semi-infinite 
% matrices written as a Toeplitz matrix plus a correction in the 
% upper leftmost corner
% Each matrix as well as the solution G, say the block A,
% is written as Toep(am,ap) +F*G' where
% am is the nonzero part of the first column of Toep(am,ap)
% ap is the nonzero part of the first row of Toep(am,ap)
% F and G are mxk and nxk slim matrices (k<<m,n)

output = fopen('jackson.txt', 'w');

% Dati two-node Jackson Motyer-Taylor


for caso=[1:10]
   if caso==1
      l1=1;l2=0;mu1=1.5;mu2=2;p=1;q=0;
   elseif caso==2
      l1=1;l2=0;mu1=2;mu2=1.5;p=1;q=0;
   elseif caso==3
      l1=0;l2=1;mu1=1.5;mu2=2;p=0;q=1;
   elseif caso==4
      l1=0;l2=1;mu1=2;mu2=1.5;p=0;q=1;
   elseif caso==5
      l1=1;l2=1;mu1=2;mu2=2;p=0.1;q=0.8;
   elseif caso==6
      l1=1;l2=1;mu1=2;mu2=2;p=0.8;q=0.1;
   elseif caso==7
      l1=1;l2=1;mu1=2;mu2=2;p=0.4;q=0.4;
   elseif caso==8
      l1=1;l2=1;mu1=10;mu2=10;p=0.5;q=0.5;
   elseif caso==9
      l1=1;l2=5;mu1=10;mu2=15;p=0.4;q=0.9;
   elseif caso==10
      l1=5;l2=1;mu1=15;mu2=10;p=0.9;q=0.4;
   end
  am = [l1+l2+mu1+mu2, (p-1)*mu1];
  ap = [l1+l2+mu1+mu2,-l1];
  bm = [-mu2*(1-q),0];
  bp = [-mu2*(1-q),-q*mu2];
  cm = [-l2,-p*mu1];
  cp = [-l2,0];
  bF=0;bG=bF;cF=bF;cG=bF;
  aG=-mu1;aF=1; 
%%%%%
  cG=-mu1;cF=1;aG=0;aF=0;
%%%%% 20/1/17

% Normalize 
f=am(1); am=am/f;ap=ap/f;bm=bm/f;bp=bp/f;cm=cm/f;cp=cp/f;
aF=aF/f;bF=bF/f;cF=cF/f;

A0 = cqt(am,ap,aF,aG);
A1 = cqt(bm,bp,bF,bG);
Am1 = cqt(cm,cp,cF,cG);
tic
  [G, R, B0] = cr(Am1,A0,A1,10);
tempo=toc;

% Residue
residuo2=norm(Am1+A0*G+A1*G^2);

fprintf (output, '%d\t%2.2e\t%2.2e\t%d\t%d\t%d\t%d\n', caso,tempo,residuo2,length(symbol(G)),size(correction(G),1),size(correction(G),2), rank(correction(G)) );
end

