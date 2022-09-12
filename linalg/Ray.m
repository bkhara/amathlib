function [v,lam,k] = Ray(A,v0)
TOL = 1e-8;
tol = 1e10;

I = eye(size(A));
v = v0./norm(v0);
lam = v0'*A*v0;

k = 0;
while(tol>TOL && k<500)
    lamprev = lam;
    
    w = (A-lam.*I)\v;
    v = w./norm(w);
    lam = v'*A*v;
    
    tol = norm(lam-lamprev);

    k = k+1;
end
end