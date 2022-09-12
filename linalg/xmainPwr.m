%% MAIN FUNCTION POWER ITER - This calls the other two power methods
v0 = ones(5,1);
for tcase = 2
%     disp(['Case ',num2str(tcase),':']);    
    switch tcase
        case 1
            A = diag([-4,2,1,1,1])+triu(rand(5,5),1);
        case 2
            A = diag([9,2,1,5,-8])+triu(rand(5,5),1);
    end
    [v,lam,k] = Pwr1(A,v0)
    [v,lam,k] = Pwr2(A,v0)
end

%% Power Iteration  Version 1
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

%% Power Iteration  Version 2
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