% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = g on boundary (fully Dirichlet conditions)
clear;
close all;

% Vector of number of elements to be used in convergence study
mVx = 7%[20 40 50 100];
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
    
    Mx = mx+2;
    My = my+2;
    m = mx*my;
    M = Mx*My;
    
    hx = Lx/(Mx-1);
    hy = Ly/(My-1);
    
    % Prepare Grid. Calculate Mat, Force and other known data
    [XP,XPb,bcAdjacency] = generateGrid(xlimits,ylimits,m,Mx,My,M);
    NNP = size(XP,1);
    
    [A,AA,b] = setupSystem(NNP,Mx,My,hx,XP);
    %     [A,AA] = setupMat(mx,my,m,hx,hy);
    %     [f] = setupForce(m,hx,hy,XP,bcAdjacency);
    %     [uBC] = setupDirichletBC(M,m,XPb);
    %     [uan] = setupAnalyticSol(m,XP);
    %     fan = A*uan;
    %     fdif = f-fan;
    
    % Solve
    u = A\b;
    divv = 0:0.05:1;
    [xq,yq] = meshgrid(divv, divv);
    vq = griddata(XP(:,1),XP(:,2),u,xq,yq);
    surf(xq,yq,vq);
    % for development purposes
    u = ones(NNP,1);
    return;
    
    % for plot purposes
    % append the boundary condition
    U = [u;uBC];
    Uan = [uan;uBC];
    
    % error calculation
    hVx(i) = hx;
    hVy(i) = hy;
    eV(i) = norm(fdif,inf);
end

% Plot solution
if(fig)
    allptx = [XP(:,1);XPb(:,1)];
    allpty = [XP(:,2);XPb(:,2)];
    xq = linspace(xlimits(1),xlimits(2),Mx);
    yq = linspace(ylimits(1),ylimits(2),My);
    [X,Y] = meshgrid(xq,yq);
    
    fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.6 0.6]);
    subplot(1,2,1);
    Z = griddata(allptx,allpty,U,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{FDM}','Location','northwest');
    bb=title(['FDM solution for', newline,' $-\nabla^2u = -1.25e^{x+0.5y}\quad in \quad \Omega=[0,1]\times[0,1],$', newline,'$u|_{\partial \Omega} = e^{x+0.5y}$'],'FontSize',16);
    set(bb,'Interpreter','latex');
    subplot(1,2,2);
    Z = griddata(allptx,allpty,Uan,X,Y);
    surf(X,Y,Z);
    xlabel('x');
    ylabel('y');
    %legend('u_{Exact}','Location','northwest');
    bb=title(['Exact solution,', newline,' $u_{exact} = e^{x+0.5y}$'],'FontSize',16);
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

% Functions involving closed form formula9s
function [u] = analyticFunction(x,y)
u = exp(x+0.5.*y);
end
function [f] = forcingFunction(x,y)
f = -1.25.*exp(x+0.5.*y);
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

function [uBC] = setupDirichletBC(M,m,XPb)
uBC = zeros(M-m,1);
for i = 1:M-m
    uBC(i) = boundaryCondition(XPb(i,1),XPb(i,2));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAJOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===========================================
% ================ GRID =====================
% ===========================================
function [XP,XPb,bcAdjacency] = generateGrid(xlimits,ylimits,m,Mx,My,M)
XP = zeros(m,4);
XPb = zeros(M-m,2);
bcAdjacency = zeros(m,2);

x = linspace(xlimits(1),xlimits(2),Mx);
y = linspace(ylimits(1),ylimits(2),My);

ki = 0;
kb = 0;

k = 0;
for j = 1:My % this order of the loops is important
    xStartInd = 1;
    jL = (My-1)/2+1;
    %     if(j > jL)
    %         xStartInd = (Mx-1)/2+1;
    %     end
    
    for i = xStartInd:Mx
        % serial order of the node
        k = k+1;
        XP(k,1) = x(i);
        XP(k,2) = y(j);
        if(i==xStartInd || i==Mx || j==1 || j==My)
            XP(k,3) = 1;
        end
        if(j == jL)
            XP(k,4) = 1;
        end
        if(j > jL)
            XP(k,4) = 2;
        end
        
        %             % identify the nodes that are adjacent to boundary nodes
        %             if(i==2)
        %                 bcAdjacency(ki,1) = -1;
        %             elseif(i==Mx-1)
        %                 bcAdjacency(ki,1) = 1;
        %             end
        %             if(j==2)
        %                 bcAdjacency(ki,2) = -1;
        %             elseif(j==My-1)
        %                 bcAdjacency(ki,2) = 1;
        %             end
    end
end
nnp = size(XP,1);
u = ones(nnp,1);
scatter3(XP(:,1),XP(:,2),u);
end

% ===========================================
% ================ MATRIX ===================
% ===========================================
function [A,AA,b] = setupSystem(NNP,Mx,My,h,XP)
coef = [-1 -1 4 -1 -1]./h^2;
A = sparse(NNP,NNP);
b = zeros(NNP,1);
for k = 1:NNP
    if(XP(k,3))
        A(k,k) = 1;
    else
        south = k-Mx;
        north = k+Mx;
        %         reduced_Mx = (Mx-1)/2+1;
        %         if(XP(k,4)==1)
        %             south = k-Mx;
        %             north = k+reduced_Mx;
        %         elseif(XP(k,4)==2)
        %             south = k-reduced_Mx;
        %             north = k+reduced_Mx;
        %         end
        A(k,south) = coef(1);
        A(k,k-1) = coef(2);
        A(k,k) = coef(3);
        A(k,k+1) = coef(4);
        A(k,north) = coef(5);
        
        b(k) = 1;
    end
end
% AA is for development purpose
if(NNP<100)
    AA = full(A);
else
    AA = 0;
end
end

function [A,AA] = setupMat(mx,my,m,hx,hy)
diag = 2*(1/hx^2 + 1/hy^2).*ones(m,1);
offdiag_far = (-1/hy^2).*ones(m-3,1);
% smaller blocks for the near diagonal
o1 = (-1/hx^2).*ones(mx-1,1);
% Fit into the bigger vector
offdiag_near = zeros(m-1,1);
for i = 1:my
    j = (i-1)*mx + 1;
    offdiag_near(j:j+(mx-2)) = o1;
    if(i<my)
        offdiag_near(j+(mx-1)) = 0;
    end
end
% Now we have the three ingredients
% diag, offdiag_near and offdiag_far
% Let's put them side by side in a matrix B
% This matrix B will be used in calling MATLAB's spdiag
% for creating A
B = zeros(m,5);
B(1:end-3,1) = offdiag_far;
B(1:end-1,2) = offdiag_near;
B(:,3) = diag;
B(2:end,4) = offdiag_near;
B(4:end,5) = offdiag_far;

% sparse diagonal position relative to the main diagonals
% this is to be used with MATLAB's spdiag
sparse_diags = [-mx -1 0 1 mx];

% Crerate the holy sparse matrix
A = sparse(m,m);
A = spdiags(B,sparse_diags,A);

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
function [f] = setupForce(m,hx,hy,XP,bcAdjacency)
f = zeros(m,1);
g = zeros(m,1);
for i = 1:m
    x = XP(i,1);
    y = XP(i,2);
    f(i) = forcingFunction(x,y);
    if(abs(bcAdjacency(i,1)))
        g(i) = g(i) + (1/hx^2)*boundaryCondition(x + hx*bcAdjacency(i,1), y);
    end
    if(abs(bcAdjacency(i,2)))
        g(i) = g(i) + (1/hy^2)*boundaryCondition(x, y + hy*bcAdjacency(i,2));
    end
end
f=f+g;
end
