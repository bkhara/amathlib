function [L,U] = luNoPivotIKJ(A)
m = size(A,1);
for i = 2:m
    for k = 1:(i-1)
        A(i,k) = A(i,k) / A(k,k);
        for j = (k+1):m
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
        end
    end
end
U = triu(A);
L = tril(A);
L = L - diag(diag(A)) + eye(m);
end