% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = g on boundary (fully Dirichlet conditions)
% Uses the 5-POINT-STENCIL
clear;
close all;

% Vector of number of elements to be used in convergence study
mVx = [24,30,42,48,96,120];
% mVx = [64];
fig = 1;

ncases = length(mVx);
mVy = mVx;
eV = zeros(ncases,1);
hVx = zeros(ncases,1);
hVy = zeros(ncases,1);

xlimits = [0,1];
ylimits = [0,1];
Lx = xlimits(2)-xlimits(1);
Ly = ylimits(2)-ylimits(1);

% loop through all the different grid spacings in mVx
for i = 1:ncases
    mx = mVx(i);
    my = mVy(i);
    
    m = mx*my;
    
    hx = Lx/(mx-1);
    hy = Ly/(my-1);
    
    % Prepare Grid. Calculate Mat, Force and other known data
    [XP] = generateGrid(xlimits,ylimits,mx,my,m);
    [A,AA] = setupMat(mx,my,m,hx,hy);
    [f] = setupForce(mx,my,m,XP);
    %[uBC] = setupDirichletBC(M,m,XPb);
    [Uan] = setupAnalyticSol(m,XP);
    fan = A*Uan;
    fdif = f-fan;
    
    % Solve
    U = A\f;
    
    % error calculation
    hVx(i) = hx;
    hVy(i) = hy;
    eV(i) = norm(fdif,inf);
end

% Plot solution
if(fig)
    allptx = [XP(:,1)];
    allpty = [XP(:,2)];
    xq = linspace(xlimits(1),xlimits(2),mx);
    yq = linspace(ylimits(1),ylimits(2),my);
    [X,Y] = meshgrid(xq,yq);
    
    fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.6 0.4]);
    subplot(1,3,1);
    Z = griddata(allptx,allpty,U,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{FDM}','Location','northwest');
    bb=title(['FDM solution for', newline,' $-\nabla^2u = 2\pi^2\sin(\pi x)\sin(\pi y)\quad in \quad \Omega=[0,1]\times[0,1],$', newline,'$u|_{\partial \Omega} = 0$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    subplot(1,3,2);
    Z = griddata(allptx,allpty,Uan,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{Exact}','Location','northwest');
    bb=title(['Exact solution,', newline,' $u_{exact} = \sin(\pi x)\sin(\pi y)$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    subplot(1,3,3);
    Z = griddata(allptx,allpty,U-Uan,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{Exact}','Location','northwest');
    bb=title(['Error, $u - u_{exact}$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    
    
    
    %     figure;
    %     subplot(1,2,1);
    %     Z = griddata(allptx,allpty,u-uan,X,Y);
    %     surf(X,Y,Z);
end

% Plot convergence in log-scale
if(ncases > 2)
    hV = hVx;
    [p,s] = polyfit(log10(hV),log10(eV),1)
    figure;
    plot(log10(hV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(h)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|f-Au_{exact}\|_{L^{\infty}}\right)$','FontSize',16);
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
u = sin(pi*x)*sin(pi*y);
end
function [f] = forcingFunction(x,y)
f = 2*pi^2*sin(pi*x)*sin(pi*y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [XP] = generateGrid(xlimits,ylimits,mx,my,m)
XP = zeros(m,2);

x = linspace(xlimits(1),xlimits(2),mx);
y = linspace(ylimits(1),ylimits(2),my);

ki = 0;
for j = 1:my % this order of the loops is important
    for i = 1:mx
        ki = ki+1;
        XP(ki,1) = x(i);
        XP(ki,2) = y(j);
    end
end
end

% ===========================================
% ================ MATRIX ===================
% ===========================================
function [A,AA] = setupMat(mx,my,m,hx,hy)
% Crerate the holy sparse matrix
A = sparse(m,m);

for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        left = I-1; right = I+1;
        down = I-mx; up = I+mx;
        if ~(j==1 || j==my || i == 1 || i == mx)
            A(I,I) = A(I,I) + 2*(1/hx^2 + 1/hy^2);
            if i>2
                A(I,left) = A(I,left) - 1/hx^2;
            end
            if i<mx-1
                A(I,right) = A(I,right) - 1/hx^2;
            end
            if j>2
                A(I,down) = A(I,down) - 1/hy^2;
            end
            if j<my-1
                A(I,up) = A(I,up) - 1/hy^2;
            end
        end
    end
end
for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        if (i==1 || i==mx || j == 1 || j == my)
            A(I,I) = 1;
        end
    end
end

% AA is for development purpose
if(m<100)
    AA = full(A);
else
    AA = 0;
end
end

% ============================================
% ================ FORCING ===================
% ============================================
function [f] = setupForce(mx,my,m,XP)
f = zeros(m,1);
for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        x = XP(I,1);
        y = XP(I,2);
        f(I) = forcingFunction(x,y);
    end
end
for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        if (j==1 || j==my || i == 1 || i == mx)
            f(I) = 0.0;
        end
    end
end
end
