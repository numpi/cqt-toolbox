function [C,D] = align_mat(A, B, direction)
switch (direction)
    case 'rows'
        m = max(size(A,1),size(B,1));
        C = zeros(m,size(A,2));
        D = zeros(m,size(B,2));
        C(1:size(A,1),1:size(A,2)) = A;
        D(1:size(B,1),1:size(B,2)) = B;
    case 'cols'
        m = max(size(A,2),size(B,2));
        C = zeros(size(A,1),m);
        D = zeros(size(B,1),m);
        C(1:size(A,1),1:size(A,2)) = A;
        D(1:size(B,1),1:size(B,2)) = B;
    case 'both'
        [C,D] = align_mat(A,B,'rows');
        [C,D] = align_mat(C,D,'cols');
    otherwise
        error('unsupported direction for the function align_mat');
end
end
