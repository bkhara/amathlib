% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = 0 on boundary (fully Dirichlet conditions)
% solved on an L-shaped domain
clear;
close all;

% Number of internal grid points
mVx = 9%[9 19 39 79 159 319];

enable_fig = 0; % enables plot
enable_sparsity_plots = 0;

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
for icases = 1:ncases
    mx = mVx(icases);
    my = mVy(icases);
    
    Mx = mx+2;
    My = my+2;
    m = mx*my;
    M = Mx*My;
    
    hx = Lx/(Mx-1);
    hy = Ly/(My-1);
    
    % Prepare Grid. Calculate Mat, Force and other known data
    assert(rem((Mx-1),2)==0, 'Mx needs to be an odd number, so that the number of segments are even'); % Need Mx to be odd for this problem
    [XP,NNP] = generateGrid(xlimits,ylimits,Mx,My);
    [A,AA,b] = setupSystem(NNP,Mx,My,hx,XP);
    
    % Solve
    u = A\b;
    unorm = norm(u,inf); % just for a check
    
    NZ = nnz(A);
    disp([mx, NZ]);
    [L,U,P] = lu(A);
    figure; spy(L);
    figure; spy(U);
    
    % [R,pmatlab] = chol(A);
    % figure; spy(R);
    
    if(enable_sparsity_plots)
        sparsityPlot(A,mx,NNP);
    end
end

% Plot solution
if(enable_fig)
    xv = linspace(xlimits(1),xlimits(2),Mx);
    yv = linspace(ylimits(1),ylimits(2),My);
    [xq,yq] = meshgrid(xv,yv);
    
    fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.4 0.6]);
    vq = griddata(XP(:,1),XP(:,2),u,xq,yq);
    surf(xq,yq,vq);
    xlabel('x');
    ylabel('y');
    bb=title(['FDM solution (Using 5-point stencil) for', newline,' $-\nabla^2u = 1$ in an L-shaped domain with $u|_{\partial \Omega} = 0$',...
        newline, 'Grid = ', num2str(mx+2),'$\times$', num2str(mx+2)], 'FontSize',16);
    set(bb,'Interpreter','latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sparsityPlot(A,mx,NNP)
fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.4 0.6]);
spy(A);
bb=title(['Sparsity pattern of A', newline, 'N = ',num2str(mx), ', size = ', num2str(NNP),'$\times$', num2str(NNP)], 'FontSize',16);
set(bb,'Interpreter','latex');
end

% ===========================================
% ================ GRID =====================
% ===========================================
function [XP,NNP] = generateGrid(xlimits,ylimits,Mx,My)
% Assume Mx = (2 * p1 + 1) and similarly My = (2 * p2 + 1)
% Then the following few lines calculates the total number of points
% in the grid using p1 and p2
p1 = (Mx-1)/2;
p2 = (My-1)/2;
NNP = Mx * (p2+1) + (p1+1) * p2;

% allocate the XP matrix
% 1st and 2nd column of XP are the x and y values of the relevant point
% If the point is on a boundary, then the value in the 3rd column is 1,
% otherwise 0
% The value in the 4th column is 1 if the point lies on the top half
% portion of the domain (above the y = 0.5 line)
XP = zeros(NNP,4);

% create one dimensional vectors for the grid points
x = linspace(xlimits(1),xlimits(2),Mx);
y = linspace(ylimits(1),ylimits(2),My);

k = 0; % k is the counter for the number of nodoes
for j = 1:My % this order of the loops is important
    xStartInd = 1;
    jL = (My-1)/2+1;
    if(j > jL)
        xStartInd = (Mx-1)/2+1;
    end
    
    for i = xStartInd:Mx
        k = k+1; % serial order of the point
        XP(k,1) = x(i);
        XP(k,2) = y(j);
        
        % boundary indicator
        if(i==xStartInd || i==Mx || j==1 || j==My || (j==jL && i<=(Mx-1)/2+1))
            XP(k,3) = 1;
        end
        
        % indicator for the points in the top half portion
        if(j == jL)
            XP(k,4) = 1;
        end
        if(j > jL)
            XP(k,4) = 2;
        end
    end
end
end

% ===========================================
% ================ MATRIX ===================
% ===========================================
function [A,AA,b] = setupSystem(NNP,Mx,My,h,XP)
coef = [-1 -1 4 -1 -1]./h^2;
A = sparse(NNP,NNP);
b = zeros(NNP,1);
for k = 1:NNP
    % check whether this point is a boundary point
    if(XP(k,3)==1)
        A(k,k) = 1;
    else
        % this is an interior node
        % so set up the north south neighbors in the stencil
        south = k-Mx;
        north = k+Mx;
        reduced_Mx = (Mx-1)/2+1;
        if(XP(k,4)==1)
            south = k-Mx;
            north = k+reduced_Mx;
        elseif(XP(k,4)==2)
            south = k-reduced_Mx;
            north = k+reduced_Mx;
        end
        
        %finally fill in the values in the matrix
        A(k,south) = coef(1);
        A(k,k-1) = coef(2);
        A(k,k) = coef(3);
        A(k,k+1) = coef(4);
        A(k,north) = coef(5);
        
        % fill in the values in the vector
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
