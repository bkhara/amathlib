function [x,rnormArr,k] = jacobi(A,b,varargin)
if(nargin>2)
    TOL = varargin{1};
else    
    TOL = 1e-8;
end
maxSz = 1500000;
rnormArr = zeros(maxSz,1);
m = size(A,1);

% Regular splitting, where M = diagonal of A
M = diag(diag(A));
N = M-A;

k = 0;
x = zeros(m,1) ;
r = b - A*x ; % residual

K1 = M\N;
K2 = M\b;

tol = norm(r);
while (tol > TOL)
    k=k+1;
    x = K1*x+K2;
    r = b-A*x;
    tol = norm(r);
    rnormArr(k) = tol;
%     if(mod(k,100000)==0)
%         disp(k)
%         disp(tol)
%     end
end

if(k<maxSz)
    rnormArr  = rnormArr(1:k,1);
end
end