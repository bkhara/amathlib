function [P,TRI,EDG,EdgeMP,BE,ElmEdges] = MeshGen(n)

npt = n^2; % total number of points

% the lower and upper bounds
lo = 0; hi = 1;

% equally spaced seeds
x = linspace(lo,hi,n);
y = linspace(lo,hi,n);

% create the grid points
[X,Y] = meshgrid(x,y);
P = zeros(npt,2);
for i = 1:npt
    P(i,:) = [X(i),Y(i)];
end

% Delaunay triangulation
TRI = delaunay(X,Y);
ntri = size(TRI,1);
tr = triangulation(TRI,P); % triangle object to use for the edge extraction
EDG = edges(tr);
nEdges = size(EDG,1);

% Create the midpoints
EdgeMP = zeros(nEdges,2);
for ied = 1:nEdges
    nd = P(EDG(ied,:),:);
    nd = nd(1,:)+nd(2,:);
    EdgeMP(ied,:) = 0.5.*nd;
end

% Store the boundary edges
nBoundaryEdges = 4*(n-1);

BE = zeros(nBoundaryEdges,1); % The boundary edges
InternalEdgeData = zeros(nEdges-nBoundaryEdges,2);
BoundaryEdgeData = zeros(nBoundaryEdges,2);

boundaryCounter = 0;
internalCounter = 0;
for i = 1:nEdges
xy = EdgeMP(i,:);
awayFromLoB = xy > lo;
awayFromHiB = xy < hi;
summ = sum(and(awayFromLoB,awayFromHiB));

if(summ < 2) %if this is a boundary edge
    boundaryCounter = boundaryCounter + 1;
    BE(boundaryCounter,1) = i;
    BoundaryEdgeData(boundaryCounter,:) = EDG(i,:);
%     BE(counter,2:3) = EDG(i,:);
%     BE(counter,4:5) = xy;
else
    internalCounter = internalCounter + 1;
    InternalEdgeData(internalCounter,:) = EDG(i,:);
end

end

% Redefine the EDG matrix
EDG = [InternalEdgeData;BoundaryEdgeData];

% Calculate the midpoints once again
EdgeMP = zeros(nEdges,2);
for ied = 1:nEdges
    nd = P(EDG(ied,:),:);
    nd = nd(1,:)+nd(2,:);
    EdgeMP(ied,:) = 0.5.*nd;
end

% Construct the ElmEdges matrix
ElmEdges = zeros(ntri,3);
local_edges = zeros(3,2);
for i = 1:ntri
    nodelist = TRI(i,:);
    nodes = sort(nodelist);
    local_edges(1,:) = [nodes(1),nodes(2)];
    local_edges(2,:) = [nodes(1),nodes(3)];
    local_edges(3,:) = [nodes(2),nodes(3)];
    
    for j = 1:3
        nd1 = local_edges(j,1);
        nd2 = local_edges(j,2);
        temp = and(EDG(:,1)==nd1,EDG(:,2)==nd2); % this vector will have 1 at the edge position
        index = find(temp==1);
        ElmEdges(i,j) = index;
    end
end

Z = exp(-X.^2-Y.^2);
trimesh(TRI,X,Y,Z);
view(0,90);
end

