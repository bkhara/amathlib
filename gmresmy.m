function [x,k,rel] = gmresmy(A,b,xi,restart)
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
    cs = zeros(n+1,1);
    sn = zeros(n+1,1);
    beta = zeros(n+1,1);
    
    Q(:,1) = r./rnorm;
    beta(1) = rnorm;
    %##########################################
    for k = 1:n
        % Arnoldi step
        [q,h] = arnoldi(A,Q,k);
        Q(:,k+1) = q; % update the last column of Q
        
        % Apply the previous Givens rotations to this new column
        % (from 1-to-k element in h
        for i = 1:k-1
            g = [cs(i),sn(i); -sn(i),cs(i)];
            h(i:i+1) = g*h(i:i+1);
            %         temp1 = h(i);
            %         temp2 = h(i+1);
            %         h(i)   =  cs(i)*temp1 + sn(i)*temp2;
            %         h(i+1) = -sn(i)*temp1 + cs(i)*temp2;
        end
        %     Get the current sin and cos values
        [cs(k),sn(k)] = GivensRotationCoeff(h(k), h(k+1));
        
        Giv = [cs(k),sn(k); -sn(k),cs(k)];
        
        % Now apply the current rotation to the current column
        %     h(k)   = cs(k)*h(k) + sn(k)*h(k+1);
        %     h(k+1) = 0.0;
        h(k:k+1) = Giv*h(k:k+1);
        h(k+1) = 0.0; % This line is redundant, but just making sure that the value is actually zero
        H(1:k+1,k) = h; % Append this new column to the end of H
        
        % Apply the current Givens rotation to the beta vector
        %     beta(k+1) = -sn(k)*beta(k);
        %     beta(k)   = cs(k)*beta(k);
        beta(k:k+1) = Giv*beta(k:k+1);
        
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

