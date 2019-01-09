function [v,lam,k] = Pwr2(A,v0)
TOL = 1e-8;
tol = 1e10;

v = v0./norm(v0);
lam = 0;
k = 0;
while(tol>TOL && k<500)
    lamprev = lam;
    
    w = A*v;
    v = w./norm(w);
    lam = v'*A*v;
    
    tol = norm(lam-lamprev);

    k = k+1;
end
end