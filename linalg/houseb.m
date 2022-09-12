function [R,b] = houseb(R,b)
% This is Householder QR method for solving linear systems
m = size(R,1);
n = size(R,2);
% W = zeros(m,n);

for k = 1:n
    x = R(k:m,k);
    
    % Creating e1
    e1 = zeros(length(x),1);
    e1(1) = 1;
    
    % Create v (i.e., the reflection vector)
    v = sign(x(1))*norm(x).*e1 + x;
    v = v./norm(v); % normalize
    % W(k:m,k) = v; % Save v as a column of W

    % Create the reflection matrix
    F = eye(m-k+1) - 2.*v*v';
    P = eye(m);
    P(k:m,k:m) = F;
    
    % Apply the reflector
    R = P*R;
    b = P*b;
    %R(k:m,k:n) = R(k:m,k:n) - 2.*W(k:m,k)*(W(k:m,k))'*R(k:m,k:n);
end
end