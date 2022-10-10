clear;
close all;

mVx = [49 99 199 399 799 999 1199 1599 3199];
ncase = length(mVx);
kV = zeros(ncase,1);
eV = zeros(ncase,1);

for icase = 1:ncase
% mesh
m = mVx(icase);
h = 1/(m+1);
x = linspace(0,1,m);

% pde data
a = 1;

% initial data
u = zeros(m,1);
U = initialC(a,x,u,m);
Unext = zeros(size(U));

% time discretization
T = 2;
k = h/2;
nts = floor(T/k);
time = linspace(0,T,nts);

% derivative operators
g1 = [1 -6 3 2]./(6*h);
g2 = [0 1 -2 1]./(h^2);
g3 = [-1 3 -3 1]./(h^3);

G1 = -a*k.*g1;
G2 = (1/2)*((a*k)^2).*g2;
G3 = -(1/6)*((a*k)^3).*g3;

G = G1+G2+G3;

% explicit time-stepping loop
for j=1:nts
    for i=2:m
        if (i == m)
            stencil = [U(m-2), U(m-1), U(m), U(2)]';
        end
        if (i == 2)
            stencil = [U(m-1), U(m), U(2), U(3)]';
        end
        if (i == 3)
            stencil = [U(m), U(2), U(3), U(4)]';
        end
        if (i>=4 && i<=m-1)
            stencil = [U(i-2), U(i-1), U(i), U(i+1)]';
        end
        Unext(i) = U(i) + G * stencil;
    end
    U = Unext;
end

uan = zeros(m,1);
for i=1:m
    uan(i)=exactSol(a,x(i),T);
end
kV(icase) = k;
eV(icase) = norm((U-uan),inf);
end

figure;
hold on;
plot(x,uan,'g','LineWidth',2)
plot(x,U,'b*');
legend('anaylitical solution','numerical solution');

% Plot convergence in log-scale
if(ncase > 2)
    [p,s] = polyfit(log10(kV),log10(eV),1)
    [p,s] = polyfit(log10(kV(5:end)),log10(eV(5:end)),1)
    figure;
    plot(log10(kV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(k)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|u-u_{exact}\|_{L^{\infty}}\right)$','FontSize',16);
    set(bb,'Interpreter','latex');
end


function [u] = exactSol(a,x,t)
u = 2*exp(-200*(x-0.5)^2)*cos(40*pi*(x-a*t));
end

function u = initialC(a,x,u,m)
for i = 1:m
    u(i) = exactSol(a,x(i),0);
end
end
