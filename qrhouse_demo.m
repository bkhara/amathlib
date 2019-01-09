function [Q,R] = qrhouse(R)
m = size(R,1);
n = size(R,2);
W = zeros(m,n);
Q = eye(m);
for j = 1:n
    j
    Q
    R    
    x = R(j:m,j);
    
    % Creating e1
    e1 = zeros(length(x),1);
    e1(1) = 1;
    
    % Create v (i.e., the reflection vector)
    v = sign(x(1))*norm(x).*e1 + x;
    v = v./norm(v); % normalize
    W(j:m,j) = v; % Save v as a column of W

    % Create the reflection matrix
%     pp=m-j+1
%     ppp=2.*v*v'
    F = eye(m-j+1) - 2.*v*v';
    P = eye(m);
    P(j:m,j:m) = F;
    
    % Apply the reflector
    R = P*R;
    Q = Q*P';
    %R(j:m,j:n) = R(j:m,j:n) - 2.*W(j:m,j)*(W(j:m,j))'*R(j:m,j:n);
end
end