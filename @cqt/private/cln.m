function y = cln(y)
%CQT Clean the relatively small components of a vector. 
%
% Y = CLN(X) returns the thresholded vector, i.e. the vector where the
%     trailing components relatively smaller than eps have been removed. 
%
% Author: Dario Bini <bini@dm.unipi.it>

  relative = true;  
  epsilon = eps;
  
  n = length(y);
  
  if relative
    nrm = norm(y,'inf')*epsilon;
  else
    nrm=epsi;
  end
  
  for i=n:-1:1
    if abs(y(i))>nrm
      break
    end
  end
  
  y = y(1:i);

