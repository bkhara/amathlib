% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = g on boundary (fully Dirichlet conditions)
% Uses the 9-point stencil for the discretization
clear;
close all;

% Vector of number of elements to be used in convergence study
mVx = [12 48 72 96];
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
    h=hx;
    % Prepare Grid. Calculate Mat, Force and other known data
    [XP,XPb,bcAdjacency] = generateGrid(xlimits,ylimits,m,Mx,My,M);
    [A,AA,f] = setupMat(mx,my,m,h,XP);
    [uBC] = setupDirichletBC(M,m,XPb);
    [uan] = setupAnalyticSol(m,XP);
    fan = A*uan;
    fdif = f-fan;
    
    % Solve
    u = A\f;
    
    R = chol(A);
    
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
    bb=title(['FDM (9point Laplacian) solution for', newline,' $-\nabla^2u = -1.25e^{x+0.5y}\quad in \quad \Omega=[0,1]\times[0,1],$', newline,'$u|_{\partial \Omega} = e^{x+0.5y}$'],'FontSize',16);
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
    %     Z = griddata(allptx,allpty,U-Uan,X,Y);
    %     surf(X,Y,Z);
end

% Plot convergence in log-scale
if(ncases > 2)
    hV = hVx;
    [p,s] = polyfit(log10(hV),log10(eV),1);
    figure;
    plot(log10(hV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(h)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|f-Au_{exact}\|_{L^{\infty}}\right)$','FontSize',16);
    set(bb,'Interpreter','latex');
    boxdim = [.2 .5 .3 .3];
    str = strcat("Slope = ",num2str(p(1)));
    annotation('textbox',boxdim,'String',str,'FitBoxToText','on');
    %grid on;
    %hAx=gca;  % avoid repetitive function calls
    % set(hAx,'xminorgrid','on','yminorgrid','on')
    
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
f = -(5/4).*exp(x+0.5.*y);
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

% ======================================================
% ================ MATRIX and VECTOR ===================
% ======================================================
function [A,AA,b] = setupMat(mx,my,m,h,XP)
% Crerate the holy sparse matrix
A = sparse(m,m);
b = zeros(m,1);
coef = [-1 -4 -1 -4 20 -4 -1 -4 -1]./(6*h^2);
for i = 1:mx
    for j = 1:my
        k = (j-1)*mx+i;
        
        x = XP(k,1);
        y = XP(k,2);
        
        bottomLeft = (k-mx) -1;
        topLeft = (k+mx) - 1;
        
        %%%%%%%%%%% A-matrix %%%%%%%%%%%%
        % The eternal diagonal element
        A(k,k) = coef(5);
        
        % North-South-East-West in the stencil (interior points of course)
        if(i > 1)
            A(k,k-1) = coef(4);
        end
        if(i < mx)
            A(k,k+1) = coef(6);
        end
        if(j > 1)
            A(k,bottomLeft + 1) = coef(2);
        end
        if(j < my)
            A(k,topLeft + 1) = coef(8);
        end
        
        % Corner points (interior of course)
        if(i>1 && j>1)
            A(k,bottomLeft) = coef(1);
        end
        if(i<mx && j>1)
            A(k,bottomLeft + 2) = coef(3);
        end
        if(i > 1 && j < my)
            A(k,topLeft) = coef(7);
        end
        if(i < mx && j < my)
            A(k,topLeft + 2) = coef(9);
        end
        
        %%%%%%%%%%% b-vector %%%%%%%%%%%%
        % Set up the HOLY b-vector (so holy that it took me 1 hour to code
        % it right
        
        % The very first thing is to add the boundary condition to the
        % b-vector
        % first take care of the North-South-East-West points around (i,j)
        if(i==1)
            b(k) = b(k) + (4/(6*h^2))*boundaryCondition(x-h,y);
        elseif(i==mx)
            b(k) = b(k) + (4/(6*h^2))*boundaryCondition(x+h,y);
        end
        if(j==1)
            b(k) = b(k) + (4/(6*h^2))*boundaryCondition(x,y-h);
        elseif(j==my)
            b(k) = b(k) + (4/(6*h^2))*boundaryCondition(x,y+h);
        end
        
        % Now take care of the corner pointsn in the stencil
        if(~(i>1 && j>1))
            b(k) = b(k) + (1/(6*h^2))*boundaryCondition(x-h,y-h);
        end
        if(~(i>1 && j<my))
            b(k) = b(k) + (1/(6*h^2))*boundaryCondition(x-h,y+h);
        end
        if(~(i<mx && j>1))
            b(k) = b(k) + (1/(6*h^2))*boundaryCondition(x+h,y-h);
        end
        if(~(i<mx && j<my))
            b(k) = b(k) + (1/(6*h^2))*boundaryCondition(x+h,y+h);
        end
        
        % Now add the actual forcing and the Laplacian correction due to
        % the 9-point stencil
        b(k) = b(k) + forcingFunction(x,y) + (h^2/12)*laplacianOfForce(x,y);
    end
end
% AA is for development purpose
if(m<100)
    AA = full(A).*(6*h^2);
else
    AA = 0;
end
end
