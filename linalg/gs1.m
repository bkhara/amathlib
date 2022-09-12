function [Q,R] = gs1(A)
% Get the sizes
[m,n] = size(A);

% Allocate the matrices
Q = zeros(m,m);
R = zeros(m,n);

% Loop on all the columns
for j = 1:n
    vj = A(:,j);
    for i = 1:(j-1) % loop on the columns till the previous column... this loop starts with the second column
        qi = Q(:,i); % all the previous orthonormal vectors
        R(i,j) = qi'*vj; % component of vj in the qi direction
        vj = vj-R(i,j)*qi; % repeatedly subtract the components
    end
    R(j,j) = norm(vj,2);
    Q(:,j) = vj./R(j,j); % normalize
end
end