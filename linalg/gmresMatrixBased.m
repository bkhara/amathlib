function [x,k,rel] = gmresMatrixBased(A,b,xi,restart)
rtol = 1e-10;
atol = 1e-16;

m = size(A,1); % A is square
n = restart;

bnorm = norm(b);
x = xi; % initial guess

for step = 1:m
    r = b-A*x;
    rnorm = norm(r);
    
    H = zeros(n+1,n);
    Q = zeros(m,n);
    beta = zeros(n+1,1);
    
    Q(:,1) = r./rnorm;
    beta(1) = rnorm;
    U = eye(n+1);
    %##########################################
    for k = 1:n
        % Arnoldi step
        [q,h] = arnoldi(A,Q,k);
        Q(:,k+1) = q; % update the last column of Q
        H(1:k+1,k) = h; % Append this new column to the end of H
        
        % Apply the previous Givens rotations to the last column
        H(:,k) = U*H(:,k);
        % Get the current sin and cos values
        [G] = MyGivensRotation(H,k,n);
        % Update the overall rotation matrix
        U = G*U;
        % Apply the current rotation to H
        H = G*H;
        % Apply the current rotation to beta
        beta(k:k+1) = G(k:k+1,k:k+1)*beta(k:k+1);
        
        errorNorm = abs(beta(k+1)); % The last term in beta is also the error
        relativeError = errorNorm / bnorm;
        
        % legal clause for freedom
        if ( relativeError <= rtol || errorNorm <= atol)
            break;
        end
    end
    % y = H(1:k,1:k) \ beta(1:k);
    y = uptrisolver(H(1:k,1:k),beta(1:k));
    x = x + Q(:,1:k)*y;
    % legal clause for freedom
    if (relativeError <= rtol || errorNorm <= atol)
        break;
    end
end
rel = relativeError;
end

