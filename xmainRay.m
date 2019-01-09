%% MAIN FUNCTION  - Rayleigh Iteration
A = diag([9,2,1,5,-8])+triu(rand(5,5),1);
for i = 1:5
    switch i
        case 1
            v0 = zeros(5,1);
        case 2
            v0 = ones(5,1);
        case 3
            v0 = randi([1 10],5,1);
        case 4
            v0 = rand(5,1);
        case 5
            v0 = 100.*ones(5,1);
    end
    disp('v0 = ');
    disp(v0);
    [v,lam,k] = Ray(A,v0)
end

%% Function - Rayleigh Iteration
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