function [Q,R] = qrmgs(Q)
% This is modified/improved Gram-Schmidt method for QR factorization
% Get the sizes
[m,n] = size(Q);
% Allocate the matrices
R = zeros(m,n);

for i = 1:n
    R(i,i) = norm(Q(:,i),2);
    Q(:,i) = Q(:,i)./R(i,i); % by this time A(:,i) is orthogonal to all the previous column (think about it)
    for j = i+1:n
        R(i,j) = Q(:,i)'*Q(:,j); % all the coefficients corresponding to the ith basis in Q
        Q(:,j) = Q(:,j) - R(i,j).*Q(:,i); % make all the subsequent columns orthogonal to the ith column
    end
end
end