function [Q,R] = mgs(A)
% Get the sizes
[m,n] = size(A);
% Allocate the matrices
Q = zeros(m,m);
R = zeros(m,n);

% Loop over all the columns
for i = 1:n
    R(i,i) = norm(A(:,i),2);
    Q(:,i) = A(:,i)./R(i,i); % by this time A(:,i) is orthogonal to all the previous column (think about it)
    for j = i+1:n
        R(i,j) = Q(:,i)'*A(:,j); % all the coefficients corresponding to the ith basis in Q
        A(:,j) = A(:,j) - R(i,j).*Q(:,i); % make all the subsequent columns orthogonal to the ith column
    end
end
end