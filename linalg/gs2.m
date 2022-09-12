% This is the classical G-S method, but using prejector matrices
function [Q,R] = gs2(A)
% Get the sizes
[m,n] = size(A);
% Allocate the matrices
Q = zeros(m,m);
R = zeros(m,n);
I = eye(m,m);
% Loop on all the columns
for j = 1:n
    R(1:j-1,j) = Q(:,1:j-1)'*A(:,j); % the first j-1 coefficients in the jth column of R
    P = I - Q(:,1:j-1)*Q(:,1:j-1)'; % projection matrix that projects the jth column to a space that is orthogonal to al the other j-1 vectors
    vj = P*A(:,j); % the component of the jth column in the new direction
    R(j,j) = norm(vj,2); % the diagonal element
    Q(:,j) = vj./R(j,j); % normalize
end
end