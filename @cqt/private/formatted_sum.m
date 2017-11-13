function C = formatted_sum(A, B)
%FORMATTED_SUM Sum two matrices of different size, padding them with zeros.

C = zeros(max(size(A, 1), size(B, 1)), ...
        max(size(A,2), size(B,2)));

C(1:size(A,1), 1:size(A,2)) = A;
C(1:size(B,1), 1:size(B,2)) = C(1:size(B,1), 1:size(B,2)) + B;

end

