clc;
clear;
n =3;
low = 0;
hi = 1;
x = linspace(low,hi,n);
y = linspace(low,hi,n);
[X,Y] = meshgrid(x,y);
tri = delaunay(X,Y);

ntri = size(tri,1);
npt = n^2;
P = zeros(npt,2);
for i = 1:npt
    P(i,:) = [X(i),Y(i)];
end
TR = triangulation(tri,P);
edg = edges(TR);

% nEdges = size(edg,1);
elementEdges = zeros(ntri,3);
local_edges = zeros(3,2);
for i = 1:ntri
    nodelist = tri(i,:);
    nodes = sort(nodelist);
    local_edges(1,:) = [nodes(1),nodes(2)];
    local_edges(2,:) = [nodes(1),nodes(3)];
    local_edges(3,:) = [nodes(2),nodes(3)];
    
    for j = 1:3
        nd1 = local_edges(j,1);
        nd2 = local_edges(j,2);
        temp = and(edg(:,1)==nd1,edg(:,2)==nd2); % this vector will have 1 at the edge position
        index = find(temp==1);
        elementEdges(i,j) = index;
    end
end


Z = exp(-X.^2-Y.^2);
trimesh(tri,X,Y,Z);
view(0,90);