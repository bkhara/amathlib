% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = g on boundary (fully Dirichlet conditions)
% Uses the 5-POINT-STENCIL
clear;
close all;

% Vector of number of elements to be used in convergence study
% mVx = [24,30,42,48,96,120];
mVx = [64];
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
    
    % Prepare known values
    [Uan] = setupAnalyticSol(m,XP);
    [Kan] = setupDiffusivity(m,XP);
    
    [A,AA] = setupMat(mx,my,m,hx,hy,Kan);
    [f] = setupForce(mx,my,m,XP);
    %[uBC] = setupDirichletBC(M,m,XPb);
    
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
    
    fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.8 0.4]);
    
    subplot(1,3,1);
    Z = griddata(allptx,allpty,U,X,Y);
    contourf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{FDM}','Location','northwest');
    bb=title(['FDM solution for', newline,' $-\nabla^2u = 0\quad in \quad \Omega=[0,1]^2,$',...
        newline,'$u(0,y) = 1, \quad u(1,y) = 0$,',...
        newline, '$\quad (n\cdot\nabla u)|_{y=0 \cup y=1}=0$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    
    if rem(mx,2)==0
        mididx = mx/2;
    else
        mididx = (mx+1)/2;
    end
    subplot(1,3,2);
    s1 = (mididx-1)*mx + 1;
    plot(xq,U(s1:s1+mx-1));
    xlabel('x'); ylabel('u');
    bb=title(['At $y=0.5$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    subplot(1,3,3);
    plot(yq,U(mididx:mx:m));
    xlabel('y'); ylabel('u');
    bb=title(['At $x=0.5$'],'FontSize',16);
    set(bb,'Interpreter','latex');
%     KZ = griddata(allptx,allpty,Kan,X,Y);
%     surf(X,Y,KZ);
end
% save('/home/khara/Biswajit/Projects/poisson_stochastic/io/left-to-right-heat-kl/u.mat','U');
save('u.mat','U');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Analytical Formulas %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions involving closed form formulas
function [u] = analyticFunction(x,y)
u = sin(pi*x)*sin(pi*y);
end
function [k] = diffusivityFunction(x,y)
% k = 1.0 + 2.*(sin(pi*x)*sin(pi*y))^2;
eta_x = 0.5;
eta_y = 0.5;
rand_tensor_list = [-0.19202,   0.046535, -0.74652,   1.9952,   0., 0. ];
kl_sum = construct_KL_sum_2D(x, y, rand_tensor_list, eta_x, eta_y);
k = exp(kl_sum);
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
function [k] = setupDiffusivity(m,XP)
k = zeros(m,1);
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);
    k(i) = diffusivityFunction(x,y);
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
function [A,AA] = setupMat(mx,my,m,hx,hy,Kan)
% Crerate the holy sparse matrix
A = sparse(m,m);

for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        left = I-1; right = I+1;
        down = I-mx; up = I+mx;
        if ~(j==1 || j==my || i == 1 || i == mx)
            
            kij = Kan(I); % k at the current node
            kleft = Kan(left); kright = Kan(right);
            kdown = Kan(down); kup = Kan(up);
            kxij = (kright - kleft)/(2*hx);
            kyij = (kup - kdown)/(2*hy);
            
            % -(\nabla k \cdot \nabla u) term (please note the -ve sign)
            A(I,left) = A(I,left) + kxij/(2*hx);
            A(I,right) = A(I,right) - kxij/(2*hx);
            
            A(I,down) = A(I,down) + kyij/(2*hy);
            A(I,up) = A(I,up) - kyij/(2*hy);
            
            % -(k \laplacian u) term (please note the -ve sign)
            A(I,I) = A(I,I) + 2*(1/hx^2 + 1/hy^2)*kij;
            A(I,left) = A(I,left) - 1/hx^2*kij;
            A(I,right) = A(I,right) - 1/hx^2*kij;
            A(I,down) = A(I,down) - 1/hy^2*kij;
            A(I,up) = A(I,up) - 1/hy^2*kij;
        end
        if j==1
            A(I,I) = A(I,I) + 1/hy;
            A(I,up) = A(I,up) - 1/hy;
        elseif j==my
            A(I,I) = A(I,I) + 1/hy;
            A(I,down) = A(I,down) - 1/hy;
        end
    end
end
for i = [1,mx]
    for j = 1:my
        I = (j-1)*mx + i;
        A(I,:) = 0;
        A(I,I) = 1;
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
% for j = 1:my
%     for i = 1:mx
%         I = (j-1)*mx + i;
%         x = XP(I,1);
%         y = XP(I,2);
%         f(I) = forcingFunction(x,y);
%     end
% end
for j = 1:my
    for i = 1:mx
        I = (j-1)*mx + i;
        if i==1
            f(I) = 1.0;
        end
        if i==my
            f(I) = 0.0;
        end
    end
end
end
