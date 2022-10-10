% This code works for the 2-D BVP
% -(u_xx + u_yy) = f
% u = 0 on boundary (fully Dirichlet conditions)
% solved on an L-shaped domain
% This code uses a 5-pt stencil for the discretisation
clear;
close all;

% Number of internal grid points
mVx = [9 19 39 79 159 319];

enable_fig = 1; % enables plot
enable_sparsity_plots = 1;

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
    h = hx;
    
    assert(rem((mx-1),2)==0, 'mx needs to be an odd number, so that the number of segments are even'); % Need mx to be odd for this problem
    
    % Generate grid
    [XP,XPb,NP,NPb] = generateGrid(xlimits,ylimits,mx,my,hx,hy);
    
    % Set up Mat A and Vec b
    % AA is only for development purpose
    [A,AA,b] = setupSystem(NP,mx,hx,XP);
    
    % Solve the linear system
    u = A\b;
    
    unorm = norm(u,inf); % just for a check
    
    % Vector of boundary conditions
    ub = zeros(NPb,1);
    
    % Append the boundary values to the solved values (for plot purposes)
    u = [u;ub];
    
    % Sparsity plots of A and R
    % Cholesky factorization of A
    Ra = chol(A);
    
    % Reverse Cuthill-Mckee permutation
    P = symrcm(A);
    % Construct the permuted matrix
    B = A(P,P);
    % Cholesky of B
    Rb = chol(B);
    
    nzA = nnz(A);
    nzRa = nnz(Ra);
    nzB = nnz(B);
    nzRb = nnz(Rb);
    
    if(enable_sparsity_plots)
        if (mx==19 || mx==39)
            sparsityPlot(A,mx,NP,'A');
            sparsityPlot(B,mx,NP,'B');
            sparsityPlot(Ra,mx,NP,'$R_A$');
            sparsityPlot(Rb,mx,NP,'$R_B$');
        end
    end
    if(icases==1)
        disp('====================== NONZERO COUNTS ====================');
    end
    disp(['mx = ',num2str(mx),', nnzA = ',num2str(nzA),', nnzR_A = ',num2str(nzRa),', nnzB = ',num2str(nzB),', nnzR_B = ',num2str(nzRb)]);
    
    % Plot solution
    if(enable_fig && mx==39)
        xv = linspace(xlimits(1),xlimits(2),Mx);
        yv = linspace(ylimits(1),ylimits(2),My);
        [xq,yq] = meshgrid(xv,yv);
        
        fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.4 0.6]);
        vq = griddata([XP(:,1);XPb(:,1)],[XP(:,2);XPb(:,2)],u,xq,yq);
        surf(xq,yq,vq);
        xlabel('x');
        ylabel('y');
        bb=title(['FDM solution (Using 5-point stencil) for', newline,' $-\nabla^2u = 1$ in an L-shaped domain with $u|_{\partial \Omega} = 0$',...
            newline, 'Grid = ', num2str(mx+2),'$\times$', num2str(mx+2)], 'FontSize',16);
        set(bb,'Interpreter','latex');
        view(0,90);
        colorbar;
    end
end
% shp = alphaShape(XP(:,1),XP(:,2));
% plot(shp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sparsityPlot(A,mx,NNP,matName)
fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.4 0.6]);
spy(A);
bb=title(['Sparsity pattern of ',matName, newline, 'N = ',num2str(mx), ', size = ', num2str(NNP),'$\times$', num2str(NNP)], 'FontSize',16);
set(bb,'Interpreter','latex');
end

% ===========================================
% ================ GRID =====================
% ===========================================
function [XP,XPb,NP,NPb] = generateGrid(xlimits,ylimits,mx,my,hx,hy)
% XP contains the positions and info related to the internal points
% NP = number of internal points
% XPb contains the positions of the boundary points
% NPb = number of boundary points

% Assume mx = (2 * p1 + 1) and similarly my = (2 * p2 + 1)
% Then the following few lines calculates the total number of internal points
% in the grid using p1 and p2
p1 = (mx-1)/2;
p2 = (my-1)/2;
NP = mx * p2 + p1 * (p2+1);

% Assume hx = hy = h
h = hx;

% Allocate the XP matrix
% 1st and 2nd column of XP are the x and y values of the relevant point
% If the point is adjacent to an x boundary, then the value in the 3rd
% column is 1, otherwise 0
% If the point is adjacent to an y boundary, then the value in the 4th
% column is 1, otherwise 0
% The value in the 5th column is 1 if the point lies on the top half
% portion of the domain (above the y = 0.5 line)
XP = zeros(NP,5);

% create one dimensional vectors for the grid points
x = linspace(xlimits(1)+h,xlimits(2)-h,mx);
y = linspace(ylimits(1)+h,ylimits(2)-h,my);

k = 0; % k is the counter for the number of nodoes
for j = 1:my % this order of the loops is important
    xStartInd = 1;
    jL = (my-1)/2;
    if(j > jL)
        xStartInd = (mx-1)/2+2;
    end
    
    for i = xStartInd:mx
        k = k+1; % serial order of the point
        XP(k,1) = x(i);
        XP(k,2) = y(j);
        
        % Boundary Indicator
        % Look for adjacency to x-parallel boundary
        if(j==1)
            XP(k,3) = -1;
        elseif(j==my || (j==jL && i<=(mx-1)/2+2))
            XP(k,3) = 1;
        end
        % Look for adjacency to y-parallel boundary
        if(i==xStartInd)
            XP(k,4) = -1;
        elseif(i==mx)
            XP(k,4) = 1;
        end
        
        % indicator for the points in the top half portion
        if(j == jL)
            XP(k,5) = 1;
        end
        if(j > jL)
            XP(k,5) = 2;
        end
    end
end

% Set up the boundary points
Mx = mx+2;
NPb = Mx + 2 * ((Mx-1)/2+1) + 2 * ((mx-1)/2) + mx;
XPb = zeros(NPb,2);

% The rest of the code is a hack of special cases
% that would work only for this L-shaped problem where the
% angle of the L-shape is aligned with the middle line from both sides

% lower x-parallel boundary
np = Mx;
xv = linspace(xlimits(1),xlimits(2),np);
yv = zeros(1,np);
istart = 1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];

% middle x-parallel boundary
np = (Mx-1)/2+1;
xv = linspace(xlimits(1),xlimits(1)+h*(Mx-1)/2,np);
yv = 0.5.*ones(1,np);
istart = iend+1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];

% upper x-parallel boundary
np = (Mx-1)/2+1;
xv = linspace(xlimits(1)+h*((mx-1)/2+1),xlimits(2),np);
yv = ones(1,np);
istart = iend+1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];

% left y-parallel boundary
np = (mx-1)/2;
xv = zeros(1,np);
yv = linspace(ylimits(1)+h,ylimits(1)+h*(np),np);
istart = iend+1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];

% middle y-parallel boundary
np = (mx-1)/2;
xv = 0.5.*ones(1,np);
yv = linspace(ylimits(1)+h*((mx-1)/2+2),ylimits(2)-h,np);
istart = iend+1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];

% right y-parallel boundary
np = mx;
xv = ones(1,np);
yv = linspace(ylimits(1)+h,ylimits(2)-h,np);
istart = iend+1;
iend = istart + np -1;
XPb(istart:iend,:) = [xv' yv'];
end

% ===========================================
% ================ MATRIX ===================
% ===========================================
function [A,AA,b] = setupSystem(NNP,mx,h,XP)
coef = [-1 -1 4 -1 -1]./h^2;
A = sparse(NNP,NNP);
b = zeros(NNP,1);
for k = 1:NNP
    % The general indices for North and South points
    south = k-mx;
    north = k+mx;
    reduced_mx = (mx-1)/2;
    % Handle the special cases for North and South point indices
    if(XP(k,5)==1) % that is, if the point is in the bottom half region
        south = k-mx;
        north = k+reduced_mx;
    elseif(XP(k,5)==2) % that is, if the point is in the top half region
        south = k-reduced_mx;
        north = k+reduced_mx;
    end
    
    % Start filling in the values in the matrix
    % Cases pertaining to x-parallel boundaries
    if(XP(k,3) == -1)
        A(k,north) = coef(5);
    elseif(XP(k,3) == 1)
        A(k,south) = coef(1);
    else
        A(k,north) = coef(5);
        A(k,south) = coef(1);
    end
    % Cases pertaining to y-parallel boundaries
    if(XP(k,4) == -1)
        A(k,k+1) = coef(4);
    elseif(XP(k,4) == 1)
        A(k,k-1) = coef(2);
    else
        A(k,k+1) = coef(4);
        A(k,k-1) = coef(2);
    end
    % Diagonal value
    A(k,k) = coef(3);
    
    % fill in the values in the vector
    b(k) = 1;
end
% AA is for development purpose
if(NNP<100)
    AA = full(A).*h^2;
else
    AA = 0;
end
end
