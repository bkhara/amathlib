function x = qrsolve(A,b)
% Matrix sizes
m = size(A,1);
n = size(A,2);

% Do a QR factorization of A and apply the same matrices to b
[A,b] = houseb(A,b);

% Remove zero rows
if(m > n)
    A = A(1:n,:);
    b = b(1:n);
end

% At this point A is triangular
x = uptrisolver(A,b);
end