function [b] = uptrisolver(A,b)
% This routine solves Ax=b where A is necessarily UPPER-TRIANGULAR
n = size(A,2);
for j = n:-1:1
    b(j) = b(j)/A(j,j);
    for k = (j-1):-1:1
        b(k) = b(k) - A(k,j)*b(j);
    end 
end
end