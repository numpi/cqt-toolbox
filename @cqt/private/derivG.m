function Gp = derivG(s,sp,advpx)
% function Gp = derivG(s,sp)
% Compute the first derivative of G= F^p by means of the Barnett factorization
% G^p = L^{-1}U, (G^p)' = L^{-1}U'-L^{-1}L'L^{-1}U
% in input: s=[s0,s1,...,s_{p-1},1], sp=[s0', s1',...s_{p-1}',s_p']
% advpx: if true, the Advanpix toobox is used for high precision computation

% By D.A. Bini, March 3, 2021

  p = length(s)-1;
  if p>1
     if advpx
        e1 = zeros(p,1,'mp');e1(1) = mp('1'); z= mp('0')*e1;
     else
        e1 = zeros(p,1);e1(1) = 1; z=0*e1;
     end
     L = toeplitz(s(p+1:-1:2),e1);
     U = toeplitz(s(1)*e1,s(1:p));
     Lp = toeplitz(sp(p+1:-1:2),z);
     Up = toeplitz(sp(1)*e1,sp(1:p));
     Li = inv(L);
     Gp = Li*Up-Li*Lp*Li*U;
  else
      if advpx
          L = [mp('1')];
          U = [s(1)];
          Lp = [mp('0')];
          Up = [sp(1)];
          Li = [mp('1')];
      else
          L = [1];
          U = [s(1)];
          Lp = [0];
          Up = [sp(1)];
          Li = [1];          
      end
      Gp = Li*Up-Li*Lp*Li*U;
  end

