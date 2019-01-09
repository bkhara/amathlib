%% MAIN FOR INVERSE ITERATION
v0 = ones(5,1);

A = diag([9,2,1,5,-8])+triu(rand(5,5),1);
mu = -8.8;
[v,lam,k] = Inv(A,v0,mu)

%% Function - Inverse Iteration
function [v,lam,k] = Inv(A,v0,mu)

TOL = 1e-8;
tol = 1e10;

I = eye(size(A));
v = v0./norm(v0);
lam = 0;
k = 0;
while(tol>TOL && k<500)
    lamprev = lam;
    
    w = (A-mu.*I)\v;
    v = w./norm(w);
    lam = v'*A*v;
    
    tol = norm(lam-lamprev);

    k = k+1;
end
end