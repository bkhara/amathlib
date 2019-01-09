function [q,h] = arnoldi(A,Q,k)
%  k should be <=(m-1)
%  Recurrence relation:
%  A * Q_{k} = Q_{k+1} * H_{k} (THUS consequently, H_{k} = Q_{k}' * A *  Q_{k})
%  Constructs the (k+1)-dimensional Krylov subspace Q_{k+1}
%  and the k-th column of the Hessenberg matrix H_{k}
%  Returns q = Q(1:m,k+1) and h = H(1:k+1,k)


q = A*Q(:,k);
h = zeros(k+1,1);
for i = 1:k % loop on rows of H
    h(i) = Q(:,i)' * q; % MGS step: calculate the coefficient in the Q(:,j) direction
    q = q - h(i) .* Q(:,i); % MGS step: remove the component that is in the Q(:,j) direction
end
qnorm = norm(q); % at this point, v is orthogonal to all the previous members of Q
if(qnorm > 1e-12)
    h(k+1) = qnorm;    
    q = q./qnorm; % normalize
else
    return;
end

end