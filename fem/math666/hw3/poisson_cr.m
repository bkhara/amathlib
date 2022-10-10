% clc;
clear;
close all;

nelV = [10]% 20 40]% 80 100];
numcases = length(nelV);
hV = 1./nelV;
l2errV = zeros(1,numcases);

for icases = 1:numcases
    
n = nelV(icases)+1; % number of (corner) nodes in one direction

% Generate the mesh data
[P,TRI,EDG,EdgeMP,BE,ElmEdges] = MeshGen(n);

% Some important data from the mesh
ntri = size(TRI,1);
nEdges = size(EDG,1);
nBoundaryEdg = size(BE,1);
NIE = nEdges - nBoundaryEdg;

% Initialize the matrix and the vector
A = sparse(NIE,NIE);
b = zeros(NIE,1);

% Assemble the global stiffness and load vector
for iel = 1:ntri
    tt = TRI(iel,1:3); % corner nodes
    cornerNodes = P(tt, 1:2); % position of corner nodes
    
    LocEdg = ElmEdges(iel,:); % Local Edges in this element
    midNodes = EdgeMP(LocEdg,:); % the midpoints of the edges
    
    [area, phi, gphi] = ElementFEMData(cornerNodes,midNodes);
    
    for i = 1:3 % loop on rows
        for j = 1:3 % loop on columns
            aa = LocEdg(i); % index of the global edge
            bb = LocEdg(j); % index of the global edge
            if (aa <=NIE && bb <=NIE)
                A(aa,bb) = A(aa,bb) + area * gphi(:,i)' * gphi(:,j);
            end
        end
        if (aa <= NIE)
            [f,~] = force(EdgeMP(aa,1:2));
            b(aa) = b(aa) + (1/3) * area * f;%gphi(:,i)' * gphi(:,j);
        end
    end
end

% Solve the linear system
u=A\b;
% Append the homogeneous boundary values
u = [u; zeros(nBoundaryEdg,1)];

% Assemble the global stiffness and load vector
errL2 = 0;
for iel = 1:ntri
    tt = TRI(iel,1:3); % corner nodes
    cornerNodes = P(tt, 1:2); % position of corner nodes
    
    LocEdg = ElmEdges(iel,:); % Local Edges in this element
    midNodes = EdgeMP(LocEdg,:); % the midpoints of the edges
    
    [area, phi, gphi,gpx,gpy] = ElementFEMData(cornerNodes,midNodes);
    
    nodal_sol = u(LocEdg);
    usol = phi'*nodal_sol;
    [~,uan] = force([gpx,gpy]);
    
    elmErr = (uan-usol)^2 * area;
    errL2 = errL2 + elmErr;
end
l2errV(icases) = sqrt(errL2);
end

% Calculate norms
logh = -log10(hV);
loge = log10(l2errV);
slope2 = -2.*logh - 1/5;
[pl2,s] = polyfit(logh, loge, 1);

% plot error norms
fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.45 0.35]);
subplot(1,1,1);
hold on;
plot(logh, loge, '*-b', 'LineWidth', 2);
plot(logh, slope2, '--b', 'LineWidth', 1);
xlabel('-log(h)');
ylabel('log(e_{L2})');
legend("||(u-u_h)||_{L2}","Slope = -2.0");

% Plot the solution
[X,Y,U] = PlotSolution(n,EdgeMP,u);

function [area, phi, gphi, gpx, gpy] = ElementFEMData(cornerNodes,midNodes)
% Construct the van der Monde matrix based on the corner nodes
M = [ones(3,1),cornerNodes]; % van der Monde matrix
% Calculate the area of the triangle using the van der Monde matrix
area = 0.5*det(M);

x1 = midNodes(1,1);
x2 = midNodes(2,1);
x3 = midNodes(3,1);

y1 = midNodes(1,2);
y2 = midNodes(2,2);
y3 = midNodes(3,2);

gpx = (x1+x2+x3)/3;
gpy = (y1+y2+y3)/3;

% Calculate the van der Monde matrix based on the edge mid points
M = [1 x1 y1;...
    1 x2 y2;...
    1 x3 y3];
Minv = inv(M);

rv = [1 gpx gpy]';
phi = Minv'*rv; % size (3x1)

% Gradients of the basis functions
gphi = Minv(2:end,:); % size (2 x 3)
end


function [f,u] = force(xy)
m = 2;
u = sin(m*pi*xy(1))*sin(m*pi*xy(2));
f = 2*m^2*pi^2*sin(m*pi*xy(1))*sin(m*pi*xy(2));
end

function [X,Y,U] = PlotSolution(n,xyv,u)
xv = xyv(:,1);
yv = xyv(:,2);

lo = 0; hi = 1;
% equally spaced seeds
x = linspace(lo,hi,n);
y = linspace(lo,hi,n);

% create the grid points
[X,Y] = meshgrid(x,y);

U = griddata(xv,yv,u,X,Y);

zeroV = zeros(n,1);
U(1,:) = zeroV';
U(n,:) = zeroV';
U(:,1) = zeroV;
U(:,n) = zeroV;

fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.45 0.35]);
subplot(1,1,1);
surf(X,Y,U);
xlabel('x');
ylabel('y');
zlabel('u');
end