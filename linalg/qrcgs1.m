function [Q,R] = qrcgs1(Q)
% This is the classical G-S method, WITHOUT using projector matrices
% Get the sizes
[m,n] = size(Q);

% Allocate the matrices
R = zeros(m,n);

% Loop on all the columns
for j = 1:n
    vj = Q(:,j);
    for i = 1:(j-1) % loop on the columns till the previous column... this loop starts with the second column
        qi = Q(:,i); % all the previous orthonormal vectors
        R(i,j) = qi'*vj; % component of vj in the qi direction
        vj = vj-R(i,j)*qi; % repeatedly subtract the components
    end
    R(j,j) = norm(vj,2);
    Q(:,j) = vj./R(j,j); % normalize
end
end