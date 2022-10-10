% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = g on boundary (fully Dirichlet conditions)
% Uses the 9-point stencil for the discretization
clear;
close all;

mVx = [20 50 80 100 150 158 200];
ncases = length(mVx);
fig=1;
mVy = mVx;
eV = zeros(ncases,1);
hVx = zeros(ncases,1);
hVy = zeros(ncases,1);

xlimits = [0,1];
ylimits = [0,1];
Lx = xlimits(2)-xlimits(1);
Ly = ylimits(2)-ylimits(1);

for i = 1:ncases
    mx = mVx(i);
    my = mVy(i);
    
    Mx = mx+2;
    My = my+2;
    m = mx*my;
    M = Mx*My;
    
    hx = Lx/(Mx-1);
    hy = Ly/(My-1);
    
    % Prepare Grid. Calculate Mat, Force and other known data
    [XP] = generateGrid(xlimits,ylimits,Mx,My,M);
    [A,b] = setupMatVec(Mx,My,M,hx,hy,XP);
    [uan] = setupAnalyticSol(M,XP);
    fan = A*uan;
    fdif = b-fan;
    
    if(ncases==1)
        h=hx;
        AA=full(A)*(6*h^2);
    end
    % Solve
    u = A\b;
    
    % error calculation
    hVx(i) = hx;
    hVy(i) = hy;
    eV(i) = norm(fdif,inf);
end

if(fig)    
    allptx = XP(:,1);
    allpty = XP(:,2);
    xq = linspace(xlimits(1),xlimits(2),mx);
    yq = linspace(ylimits(1),ylimits(2),my);
    [X,Y] = meshgrid(xq,yq);
    
    figure;
    subplot(1,2,1);
    Z = griddata(allptx,allpty,u,X,Y);
    surf(X,Y,Z);
    subplot(1,2,2);
    Z = griddata(allptx,allpty,uan,X,Y);
    surf(X,Y,Z);
    
%     figure;
%     subplot(1,2,1);
%     Z = griddata(allptx,allpty,u-uan,X,Y);
%     surf(X,Y,Z);
end

if(ncases>2)
    hV = hVx;
    [p,s] = polyfit(log10(hV),log10(eV),1)
    figure;
    plot(log10(hV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(h)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|f-A*u_{exact}\|\right)$','FontSize',16);
    set(bb,'Interpreter','latex');
end

% Debugging related things
% fan = AA*uan;
% hyb = [f fan]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Analytical Formulas %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions involving closed form formulas
function [u] = analyticFunction(x,y)
u = exp(x+0.5.*y);
end
function [f] = forcingFunction(x,y)
f = -1.25.*exp(x+0.5.*y);
end
function [d2f] = laplacianOfForce(x,y)
d2f = -(25/16).*exp(x+0.5.*y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u] = boundaryCondition(x,y)
u = analyticFunction(x,y);
end

function [u] = setupAnalyticSol(m,XP)
u = zeros(m,1);
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);
    u(i) = analyticFunction(x,y);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAJOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===========================================
% ================ GRID =====================
% ===========================================
function [XP] = generateGrid(xlimits,ylimits,Mx,My,M)
XP = zeros(M,2);

x = linspace(xlimits(1),xlimits(2),Mx);
y = linspace(ylimits(1),ylimits(2),My);

for j = 1:My % this order of the loops is important
    for i = 1:Mx
        k = (j-1)*Mx + i;        
        % serial order of the node
        % Set up the position matrices
        XP(k,1) = x(i);
        XP(k,2) = y(j);
    end
end
end

% ===========================================
% ====== A,b SETUP IN THE FEM STYLE======
% ===========================================
function [A,b] = setupMatVec(Mx,My,M,hx,hy,XP)
h=hx;
coef = [-1 -4 -1 -4 20 -4 -1 -4 -1]./(6*h^2);
A = sparse(M,M);
b = zeros(M,1);
for i = 1:Mx
    for j = 1:My
        k = (j-1)*Mx+i;
        
        x = XP(k,1);
        y = XP(k,2);
        
        if(i==1 || i==Mx || j==1 || j==My)
            A(k,k) = 1;
            b(k) = boundaryCondition(x,y);
        else
            bottomLeft = (k-Mx) -1;
            topLeft = (k+Mx) - 1;
            
            A(k,bottomLeft) = coef(1);
            A(k,bottomLeft + 1) = coef(2);
            A(k,bottomLeft + 2) = coef(3);
            
            A(k,k-1) = coef(4);
            A(k,k) = coef(5);
            A(k,k+1) = coef(6);
            
            A(k,topLeft) = coef(7);
            A(k,topLeft + 1) = coef(8);
            A(k,topLeft + 2) = coef(9);
            
            b(k) = forcingFunction(x,y) + (h^2/12)*laplacianOfForce(x,y);
        end
    end
end
end