function [x,rnormArr,k] = cg(A,b,varargin)
if(nargin>2)
    TOL = varargin{1};
else    
    TOL = 1e-8;
end
maxSz = 500;
rnormArr = zeros(maxSz,1);
m = size(A,1);

k = 0;
x = zeros(m,1) ;
r = b - A*x ; % residual
p = r; % search direction
tol = norm(r);
while (tol > TOL)
    k=k+1;
    
    Ap = A*p;
    a0 = (r' * r) / (p' * Ap) ;
    x = x + a0 .* p;
    
    r0 = r; % save the old r
    r = r - a0.*Ap; % new residual
    b0 = (r'*r)/(r0'*r0); % improvement this step
    p = r + b0.*p; % search direction improvement
    
    tol = norm(r);
    rnormArr(k) = tol;
end
if(k<maxSz)
    rnormArr  = rnormArr(1:k,1);
end
end

% function [x1] = cg(A,b)
% TOL = 1e-6;
% m = size(A,1);
% k = 0;
% x0 = zeros (m, 1 ) ;
% r0 = b - A* x0 ; % residual
% while (norm(r0) > TOL)
%     k=k+1;
%     Ar0 = A*r0;
%     a0 = (r0' * Ar0) / (Ar0' * Ar0) ;
%     x1 = x0 + a0 * r0;
%     x0 = x1;
%     r0 = b - A*x0;
% %     Orth1_Residual(k)=norm(r0) ;
% end
% end