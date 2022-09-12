function [L,U] = ilubasic0(A)
m = size(A,1);


for k = 1:(m-1)
    for i = (k+1):m
        if(A(i,k)>1e-10)
            A(i,k) = A(i,k)./A(k,k);
        end
        for j = k+1:m
            if(A(i,j)>1e-10)
                A(i,j) = A(i,j) - A(i,k).*A(k,j);
            end
        end
    end
end
L = eye(m);
U = triu(A);
L = tril(A);
L = L - diag(diag(A)) + eye(m);
end