function w = range(A)
% RANGE Plot the range of the symbol of a CQT matrix. 
%
% W = RANGE(A) computes and plot W = a(T), where T is the unit circle, and 
% a(z) the symbol of the CQT matrix A. 

% By D.A. Bini, March 3, 2021

  [am, ap] = symbol(A);

  smp = 2^20; 
  nm = length(am); np = length(ap);
  n = max(smp,2*(nm+np));
  
% check size
  if size(am,2)>1
     am = am.';
  end
  if size(ap,2)>1
     ap = ap.';
  end
  
% sample the function values
  f = zeros(n,1);
  f(1:np) = ap; 
  f(end:-1:end-nm+2) = am(2:nm);
  w = fft(f);
% plot
  plot(w,'.','linewidth',0.5,'color',[0.6,0.8,1]);

