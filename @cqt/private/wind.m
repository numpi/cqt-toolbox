function w = wind(am,ap,x,advpx)
% function w = wind(am,ap,x,advpx)
% w is the winding number of the Laurent polynomial a(z)-x
% where a(z) has coefficients in am and ap
% the winding number is evaluated by counting the roots 
% of modulus less than 1 by means of Graeffe's method
% advpx: if true, the toolbox advanpix is used for high precision computation
% By D.A. Bini, January, 2022

  gsteps = 18;
  am(1) = am(1)-x;
  ap(1) = ap(1)-x;
  m = length(am);
  n = length(ap);
  d = m+n-1;
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
  [~,vind] = max(abs(a)); ind = min(vind);
  nz = length(a)-ind;
  m = length(am)-1;
  w = nz-m;
  fprintf('Warning: in wind, max number of Graeffe steps reached\n');
  return
end

