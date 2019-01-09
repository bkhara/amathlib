function [U,H] = hessenberg(A)
% Transforms a general matrix to upper Hessenberg form by a series of 
% similarity transformation
% The orthogonalization is done by Householder reflectors
% This is a direct method, as in the number of loops are known
% Returns U and H such that
% A = UHU'
% OR
% H = U'AU
m = size(A,1);
H=A;
U=eye(size(H));
for k = 1:m-2
    x = H(k+1:m,k);
    
    % Creating e1
    e1 = zeros(length(x),1);
    e1(1) = 1;
    
    % Create v (i.e., the reflection vector)
    v = sign(x(1))*norm(x).*e1 + x;
    v = v./norm(v); % normalize

    % Create the reflection matrix
    F = eye(m-k) - 2.*v*v';
    P = eye(m);
    P(k+1:m,k+1:m) = F;
    
    % Apply the reflector by pre- and post-multiplying
    H = P*H*P';
    
    % Update (or construct) U
    U = P*U;
end
% Comment:
% U is the matrix that operates on A to create H, 
% thus it makes sense to output U
% But due to the common definiton of H given as A = U*H*U',
% here, the U is being transposed
U = U';
end