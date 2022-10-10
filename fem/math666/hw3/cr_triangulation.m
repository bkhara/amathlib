clc;
clear;

low = 0;
hi = 1;

nsegments1 = 3;
n1 = nsegments1 + 1;

x = linspace(low,hi,n1);
y = linspace(low,hi,n1);
[X1,Y1] = meshgrid(x,y);
tri = delaunay(X1,Y1);

% Z = exp(-X1.^2-Y1.^2);
% trimesh(tri,X1,Y1,Z);

nsegments2 = 2*nsegments1;
n2 = nsegments2 + 1;
x = linspace(low,hi,n2);
y = linspace(low,hi,n2);
[X2,Y2] = meshgrid(x,y);

tri2 = zeros(size(tri));
for iel = 1:size(tri,1)
    nodelist = tri(iel,:);
    xp = X1(nodelist);
    yp = Y1(nodelist);
    
    %orient = orientation(xp,yp);
    [tri_orient, new_order] = orientation(X,Y,nodelist);
    if tri_orient % upright triangle
        corner_node = new_order(1);
        row_id = rem(corner_node,n1);%corner_node - n1*column_id;
        column_id = (corner_node - row_id) / n1;
        
        
        row_on_finer_grid = 
    elseif ~tri_orient % downright triangle
        
    end
end
