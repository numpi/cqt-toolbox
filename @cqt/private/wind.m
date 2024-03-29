function w = wind(am,ap,x,advpx)
% function w = wind(am,ap,x,advpx)
% w is the winding number of the Laurent polynomial a(z)-x
% where a(z) has coefficients in am and ap
% the winding number is evaluated by counting the roots 
% of modulus less than 1 by means of Graeffe's method
% advpx: if true, the toolbox advanpix is used for high precision computation
% By D.A. Bini, March 30, 2022

  gsteps = 18;  verbose = false; %verbose = true;
  am(1) = am(1)-x;
  ap(1) = ap(1)-x;
  m = length(am);
  n = length(ap);
  d = m+n-1;
  if advpx
      epsi = mp(10^(-mp.Digits+3));
  else
      epsi = eps*10^3;    
  end
% preliminaries
  s = (-1).^(1:d)';
  s = s*sign(s(m));
  a=[ap(end:-1:1);am(2:end)];
% normalize a
  a = a/max(abs(a));
% check number of roots before starting Graeffe
  if sum(abs(a))<2
     [~,ind] = max(abs(a));
     nz = length(a)-ind;
     m = length(am)-1;
     w = nz-m;        
     return
  end
% start Graeffe 
  for i=1:gsteps 
     t = conv(a,s.*a);
     [~,ind] = max(abs(t));
     t = t/t(ind);
     a = t(1:2:length(t));
     if length(a)<d
         if advpx
              a=[zeros(d-length(a),1,'mp');a];
         else
              a=[zeros(d-length(a),1);a];
         end
     end        
         if sum(abs(a))<2
              [~,ind] = max(abs(a));
              nz = length(a)-ind;
              m = length(am)-1;
              w = nz-m;           
              return
         end
  end
  if verbose
      fprintf('Warning: in wind, max number of Graeffe steps reached\n');
  end
  
% filtering coefficients to avoid malfunctioning of roots 
  b=abs(a)>1.e-60;
  a = a.*b;
  absr = abs(roots(a));
  if advpx
      chk = mp('1')-absr;
  else
      chk =1-absr;
  end
  nz = sum(chk>epsi*2^18); % numerical correction 
  m = length(am)-1;
  w = nz-m;   
  return
end

