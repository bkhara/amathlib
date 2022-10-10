% This code works for the heat equation
% u_t = (u_xx + u_yy)
% u = g on boundary (fully Dirichlet conditions)
% Uses the 5-POINT-STENCIL
clear;
close all;

% Vector of number of elements to be used in convergence study
mVx = [5 10 20 50 100 200];%[5 10 20 50 100];%[24,30,42,48,96,120,150,200,300];
nVt = mVx;
fig = 1;

ncases = length(nVt);
mVy = mVx;
eV = zeros(ncases,1);
hVx = zeros(ncases,1);
hVy = zeros(ncases,1);
kV = zeros(ncases,1);

xlimits = [0,1];
ylimits = [0,1];
tlimits = [0,1e-3];
Lx = xlimits(2)-xlimits(1);
Ly = ylimits(2)-ylimits(1);
T = tlimits(2)-tlimits(1);

% loop through all the different grid spacings in mVx
for i = 1:ncases
    mx = mVx(i);
    my = mVy(i);
    nts = nVt(i);
    
    disp(['Nx = ',num2str(mx)]);
    
    Mx = mx+2;
    My = my+2;
    m = mx*my;
    M = Mx*My;
    
    hx = Lx/(Mx-1);
    hy = Ly/(My-1);
    h=hx;
    
    k = T / nts;
    
    % Prepare Grid. Calculate Mat, Force and other known data
    [XP,XPb,bcAdjacency] = generateGrid(xlimits,ylimits,m,Mx,My,M);
    [Axleft,Axright,Ayleft,Ayright,Aleft,Aright] = setupMat(mx,my,m,hx,hy,k);
    u = zeros(m,1);
    u = initialCondition(XP,m,u);
    t = 0;
    for its = 1:nts
        rhs_stage1 = Ayright*u;
        rhs_stage1 = setupForceTransient(t,rhs_stage1,m,hx,hy,k,XP,bcAdjacency,1);
        if(its==1)
            tic
        end
        u = Axleft \ rhs_stage1;
        
        rhs_stage2 = Axright*u;
        rhs_stage2 = setupForceTransient(t,rhs_stage2,m,hx,hy,k,XP,bcAdjacency,2);
        u = Ayleft \ rhs_stage2;
        if(its==1)
            toc
        end
        t = t + k;
    end
    
    [uBC] = setupDirichletBCTransient(M,m,XPb,T);
    [uan] = setupAnalyticSolTransient(m,XP,T);
    
    % error calculation
    hVx(i) = hx;
    hVy(i) = hy;
    kV(i) = k;
    eV(i) = norm((rhs_stage2-Ayleft*uan),inf);
    % for plot purposes
    % append the boundary condition
    u = [u;uBC];
    uan = [uan;uBC];
    
    % Plot solution
if(fig == 1 && mx==100)
    allptx = [XP(:,1);XPb(:,1)];
    allpty = [XP(:,2);XPb(:,2)];
    xq = linspace(xlimits(1),xlimits(2),Mx);
    yq = linspace(ylimits(1),ylimits(2),My);
    [X,Y] = meshgrid(xq,yq);
    
    fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.6 0.6]);
    subplot(1,2,1);
    Z = griddata(allptx,allpty,u,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{FDM}','Location','northwest');
    bb=title(['FDM solution'], 'FontSize',16);
    set(bb,'Interpreter','latex');
    subplot(1,2,2);
    Z = griddata(allptx,allpty,uan,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{Exact}','Location','northwest');
    bb=title(['Exact solution'],'FontSize',16);
    set(bb,'Interpreter','latex');
end
end



% Plot convergence in log-scale
if(ncases > 2)
    hV = hVx;
    %     [p,s] = polyfit(log10(hV),log10(eV),1)
    %     figure;
    %     plot(log10(hV),log10(eV),'bo-')
    [p,s] = polyfit(log10(kV),log10(eV),1)
    figure;
    plot(log10(kV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(h)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|f-Au_{exact}\|_{L^{\infty}}\right)$','FontSize',16);
    set(bb,'Interpreter','latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Analytical Formulas %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = initialCondition(XP,m,u)
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);
    u(i) = analyticFunctionTransient(x,y,0);
end
end

function [u] = analyticFunctionTransient(x,y,t)
m=3;
u = exp(-32*pi^2*t).*cos(4.*pi.*x).*cos(4.*pi.*y);
% u = exp(-2*(m*pi)^2*t).*sin(m*pi.*x).*sin(m*pi.*y);
% u = exp(-2*(m*pi)^2*t).*cos(m*pi.*x).*cos(m*pi.*y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u] = boundaryConditionTransient(x,y,t)
u = analyticFunctionTransient(x,y,t);
end

function [u] = setupAnalyticSolTransient(m,XP,t)
u = zeros(m,1);
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);
    u(i) = analyticFunctionTransient(x,y,t);
end
end

function [uBC] = setupDirichletBCTransient(M,m,XPb,t)
uBC = zeros(M-m,1);
for i = 1:M-m
    uBC(i) = boundaryConditionTransient(XPb(i,1),XPb(i,2),t);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAJOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===========================================
% ================ GRID =====================
% ===========================================
function [XP,XPb,bcAdjacency] = generateGrid(xlimits,ylimits,m,Mx,My,M)
XP = zeros(m,2);
XPb = zeros(M-m,2);
bcAdjacency = zeros(m,2);

x = linspace(xlimits(1),xlimits(2),Mx);
y = linspace(ylimits(1),ylimits(2),My);

ki = 0;
kb = 0;
for j = 1:My % this order of the loops is important
    for i = 1:Mx
        
        % serial order of the node
        % k = (j-1)*My + i;
        if(i==1 || i==Mx || j==1 || j==My)
            kb = kb+1;
            XPb(kb,1) = x(i);
            XPb(kb,2) = y(j);
        else
            % Set up the position matrices
            ki = ki+1;
            XP(ki,1) = x(i);
            XP(ki,2) = y(j);
            
            % identify the nodes that are adjacent to boundary nodes
            if(i==2)
                bcAdjacency(ki,1) = -1;
            elseif(i==Mx-1)
                bcAdjacency(ki,1) = 1;
            end
            if(j==2)
                bcAdjacency(ki,2) = -1;
            elseif(j==My-1)
                bcAdjacency(ki,2) = 1;
            end
        end
    end
end
end

% ===========================================
% ================ MATRIX ===================
% ===========================================
function [Axleft,Axright,Ayleft,Ayright,Aleft,Aright] = setupMat(mx,my,m,hx,hy,k)
idiag = ones(m,1);
diagx = -2*(1/hx^2).*(k/2).*ones(m,1);
diagy = -2*(1/hy^2).*(k/2).*ones(m,1);
offdiag_far = (1/hy^2).*(k/2).*ones(m-3,1);
% smaller blocks for the near diagonal
o1 = (1/hx^2).*(k/2).*ones(mx-1,1);
% Fit into the bigger vector
offdiag_near = zeros(m-1,1);
for i = 1:my
    j = (i-1)*mx + 1;
    offdiag_near(j:j+(mx-2)) = o1;
    if(i<my)
        offdiag_near(j+(mx-1)) = 0;
    end
end

Bx = zeros(m,3);
Bx(1:end-1,1) = offdiag_near;
Bx(2:end,3) = offdiag_near;

By = zeros(m,3);
By(1:end-3,1) = offdiag_far;
By(4:end,3) = offdiag_far;
% sparse diagonal position relative to the main diagonals
% this is to be used with MATLAB's spdiag
sparse_diagsx = [-1 0 1];
sparse_diagsy = [-mx 0 mx];

Axleft = sparse(m,m);
Axright = sparse(m,m);
Ayleft = sparse(m,m);
Ayright = sparse(m,m);

Bx(:,2) = diagx - idiag;
Axleft = spdiags(Bx,sparse_diagsx,Axleft);
Axleft = (-1).*Axleft;
Bx(:,2) = diagx + idiag;
Axright = spdiags(Bx,sparse_diagsx,Axright);

By(:,2) = diagy - idiag;
Ayleft = spdiags(By,sparse_diagsy,Ayleft);
Ayleft = (-1).*Ayleft;
By(:,2) = diagy + idiag;
Ayright = spdiags(By,sparse_diagsy,Ayright);

% AA is for development purpose
if(m<100)
    Aleft = full(Axleft);
    Aright = full(Axright);
else
    Aleft = 0;
    Aright = 0;
end
end

% ============================================
% ================ FORCING ===================
% ============================================
function [f] = setupForceTransient(t,f,m,hx,hy,k,XP,bcAdjacency,stage)
g = zeros(m,1);
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);    
    if(abs(bcAdjacency(i,1)))
        if(stage==1)
            g(i) = g(i) + (k/2).*(1/hx^2)*boundaryConditionTransient(x + hx*bcAdjacency(i,1), y,(t+k/2));
        elseif(stage==2)
            g(i) = g(i) + (k/2).*(1/hx^2)*boundaryConditionTransient(x + hx*bcAdjacency(i,1), y,(t+k/2));
        end
    end
    if(abs(bcAdjacency(i,2)))
        if(stage==1)
            g(i) = g(i) + (k/2).*(1/hy^2)*boundaryConditionTransient(x, y + hy*bcAdjacency(i,2),(t));
        elseif(stage==2)
            g(i) = g(i) + (k/2).*(1/hy^2)*boundaryConditionTransient(x, y + hy*bcAdjacency(i,2),(t+k));
        end
    end
end
f=f+g;
end
