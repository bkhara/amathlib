function [v,lam,k] = Pwr1(A,v0)
TOL = 1e-8;
tol = 1e10;

v = v0./norm(v0);
lam = 0;
k = 0;
while(tol>TOL && k<500)
    vprev = v;
    
    w = A*v;
    v = w./norm(w);
    lam = v'*A*v;
    
    tol = norm(v-vprev);

    k = k+1;
end
end