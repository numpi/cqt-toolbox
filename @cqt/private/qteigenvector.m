function v = qteigenvector(n,V,y,m,algo,advpx)
% v = function qteigenvector(n,V,y,m,algo,advpx)
% Compute n block components of the eigenvector v of a QT matrix
% given the parameters V, y, m, algo provided by the function qteig
% V : either Vandermonde pxp or G according to the algorithm
% y : such that WVy=0
% m : length(am)-1
% algo: algorithm
% advpx: if true, the toolbox Advanpix is used for high precision arithmetic

% By D.A. Bini, Jan, 2022

if size(y,2)>1 
   y = y.';
end
p = length(y); 
if advpx
   v = zeros(n*p,1,'mp');
else
   v = zeros(n*p,1);
end
if mod(algo,2)==1   % Vandermonde
   r = V(2,:).';
   z = r.^m;         
   z = y.*z;   
   v(1) = sum(z);               
   for i=2:n*p
      z = z.*r;
      v(i) = sum(z);
   end
else                % Frobenius
   g = V(p+1,:).'; G = V(p+1:2*p,:);
   v(1:p) = y; q=ceil(m/p); n=n+q; %%%%%%%%%%%%%%%%%%%%%%%
   for i=1:n
      v(i*p+1:i*p+p) = G*v((i-1)*p+1:i*p);
   end
   v = v(m+1:min(m+n*p,n*p+p));
end
v = v/norm(v);
end