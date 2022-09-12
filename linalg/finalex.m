%% MATH562/FINAL EXAM/Biswajit Khara
clear;
m=500; % Number of internal grid points
tol = 1e-10;
maxit = 100;
% Domain distributed force
f = 1;
% Set up the data related to the grid
[X,h] = SetUpGridData(m);
% Calculate the A and b matrices
ep = 0;
[A,b] = SetUpMatricesAB(m,h,f,ep);
% Solve the system
% disp('Matlab Cholesky');
tic
[um] = A\b; % using matlab mldivide
tm = toc;

tic
[ujac,rnormjac,kjac] = jacobi(A,b,tol);
tjac = toc;
tic
[ugs,rnormgs,kgs] = gseidel(A,b,tol);
tgs = toc;
tic
[ucg,rnormcg,kcg] = cg(A,b,1e-16);
tcg = toc;

PlotNorms('Jacobi',X,ujac,kjac,rnormjac,tjac,'b');
PlotNorms('Gauss-Seidel',X,ugs,kgs,rnormgs,tgs,'k');
PlotNorms('CG',X,ucg,kcg,rnormcg,tcg,'r');
fileNm = strcat('dat_',num2str(m));
save(fileNm);

%% Plot the norms against iteration number
function PlotNorms(method,X,u,k,rnorm,t,clr)
fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 1/2 1/3]);
hold on;
kk=log(10);
subplot(1,2,1);
plot(log(rnorm)./kk,'Color',clr,'LineWidth',2,'DisplayName',strcat('Convergence:'," ",method));
xlabel('No. of iter');
ylabel('log_{10} || r ||');
lgd = legend;
subplot(1,2,2);
[u] = AddBC(u);
plot(X,u,'Color',clr,'LineWidth',2,'DisplayName',strcat('Solution:'," ",method));
xlabel('x');
ylabel('u');
lgd=legend;
disp([method,':   iter = ',num2str(k),',  normR = ',num2str(rnorm(k)),',  time = ',num2str(t)]);
end

%% Set up the data related to the grid
function [X,h] = SetUpGridData(m)
xlo = 0;
xhi = 1;
gridpoints = m+2;
X = linspace(xlo,xhi,gridpoints);
% Steplength
h = (xhi-xlo)/(m+1);
end
%% Set up the A and b matrices
function [A,b] = SetUpMatricesAB(m,h,f,ep)
% Initialize matrices
A = zeros(m,m);
b = zeros(m,1);
% Assemble A
v = -[1 -2-ep 1]; % negative because of the negative sign occuring in the differential equation
A(1,1:2) = v(2:3);
A(m,m-1:m) = v(1:2);
for i = 2:m-1
    A(i,i-1:i+1) = v;
end
A = A./(h.^2);
% Assemble b
b(:) = f.* ones(m,1);
end
%% Incorporate the boundary condition into the solution vector
function [u] = AddBC(u)
u1 = 0;
u2 = 0;
u = [u1; u; u2];
end
%% Jacobi iteration code
function [x,rnormArr,k] = jacobi(A,b,varargin)
if(nargin>2)
    TOL = varargin{1};
else
    TOL = 1e-8;
end
maxSz = 1500000;
rnormArr = zeros(maxSz,1);
m = size(A,1);

D = diag(diag(A));
M = D;
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
%% Gauss-Seidel code
function [x,rnormArr,k] = gseidel(A,b,varargin)
if(nargin>2)
    TOL = varargin{1};
else
    TOL = 1e-8;
end
maxSz = 1500000;
rnormArr = zeros(maxSz,1);
m = size(A,1);

D = diag(diag(A));
E = -tril(A,-1);
M = D-E;
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
end
if(k<maxSz)
    rnormArr  = rnormArr(1:k,1);
end
end

%% Conjugate Gradient method code
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
%     if(k>=m)
%         break;
%     end
end
if(k<maxSz)
    rnormArr  = rnormArr(1:k,1);
end
end

