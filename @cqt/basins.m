function [AA,v,iter] = basins(A,n,x0,x1,y0,y1,algo,plotfig,advpx,digits)
% function [A,v,iter]  = qtbasins(am,ap,E,n,x0,x1,y0,y1,algo,plotfig,advpx,digits)
% If plotfig is true then the function draws the basins of attraction of the 
% fixed point iteration given by algorithm algo, where algo=1,2,3,4,
% applied to the matrix CQT(am,ap,E) in the range x0<re(x)<x1, y0<im(x)<y1
% On output: A is the nxn matrix of unint n: A(k,j) contains the integer s 
% such that the fixed point sequence obtained by starting from the point 
% z0 = x0+(k-1)(x1-x0)/(n-1) + i(y0+(j-1)(y1-y0)/(n-1)), k,j=1:n
% converges to the (s+1)-st eigenvalue. Moreover,
% A(k,j) = 256-1 if z0 is in a continuous set
%          255-1 if the sequence exits the component containing z0
%          254-1 if the sequence diverges or if exceeded the max number of iterations
%
% v:    vector with the isolated eigenvalues
% iter: overall number of iterations
% a picture of the basins is plotted with the following colors:
% x eig: color associated with the eig      (flag 1 or 3)
% 255 continuous set:       green           (flag 2)
% 254 out of the component: light gray      (flag 5)
% 253 p<q not an eig :      black           (flag 4)
% 253 Exceeded iter:        black           (flag 6)
% 253 Diverging seq.        black           (flag 7)

% By D.A. Bini, May 3, 2021

  [am, ap] = symbol(A);
  E = correction(A);

  verbose = true; maxit = 30; epsi=1.e-11; residual=false;
  AA = uint8(zeros(n));
  v = []; 
  
  % adjust the input size
  if size(am,2)~=1
     am = am.';
  end
  if size(ap,2)~=1
     ap = ap.';
  end

% if x0=x1 set default values for x0 and x1  
  if x0==x1
     a = [ap(end:-1:1);am(2:end)];
     mxval = sum(abs(a)) + norm(E,inf);
     x0 = -mxval*1.001; x1 = mxval; y0 = x0; y1 = x1;
  end
  
  % compute the eigenvalues 
  v = eig_all(A,  algo, 20, 1.e-12, 4, false, false, residual, advpx, digits);
  nz = length(v);
  
  % choose the color palette
   MP = rand(256,3);
   MP(256,:) = [0.4 0.9 0.4];   %  continuous set
   MP(255,:) = [0.8,0.8,0.8];   % out of the component
   MP(254,:) = [0 0 0];         % non-converging sequence

% construct the basins  
  iter = 0;
  for k=1:n
    for j=1:n
       z0 = x0+(j-1)*(x1-x0)/(n-1) + 1i*(y0+(k-1)*(y1-y0)/(n-1));
       [z, flg, ~, ~,~,it, ~, ~, ~] = eig(A,z0,algo,maxit,epsi,false,residual,advpx,digits);
       iter = iter+it;
       if verbose && j==1 &&  mod(k,10)==1
          fprintf('[k,j,flg]=[%d,%d,%d] \n',k,j,flg);
       end 
       if flg==1  || flg==3      %  isolated eigenvalue
          for l=1:nz
             if abs(z-v(l))<1.e-8 
                AA(k,j) = min(l-1,nz-1); 
                break
             end
          end
       elseif flg==2    % p>q: continuous eigenvalue
          AA(k,j) = 256-1;
       elseif flg==5    % change of component
          AA(k,j) = 255-1;
       elseif flg>=4    % non-converging sequence 
          AA(k,j) = 254-1;
       end
    end
  end  
  if plotfig
     figure; image([x0,x1],[y0,y1],AA); %title(' Basins ');
     colormap(MP);
     ax = gca;
     ax.YDir = 'normal';
  end
end
