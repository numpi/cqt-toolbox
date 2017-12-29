function A = merge_corrections(A)
%MERGE_CORRECTIONS

m = size(A, 1);
n = size(A, 2);

if size(A.U, 1) + size(A.W, 1) > n && ...
        size(A.V, 1) + size(A.Z, 1) > n
    A.U = [ [ A.U ; zeros(m-size(A.U,1),size(A.U,2)) ] , ...
           [ zeros(m-size(A.W,1),size(A.W,2)) ; A.W(end:-1:1,end:-1:1) ] ];
    A.V = [ [ A.V ; zeros(n-size(A.V,1),size(A.V,2)) ] , ...
          [ zeros(n-size(A.Z,1),size(A.Z,2)) ; A.Z(end:-1:1,end:-1:1) ] ];
    A.W = [];
    A.Z = [];
end

