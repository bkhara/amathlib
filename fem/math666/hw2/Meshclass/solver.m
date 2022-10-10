clear;
close all;

icase = 1;

switch icase
    case 0
        load ('squareholemesh.mat', 'p','t','NIN');
    case 1
        load ('starmesh.mat', 'p','t','NIN');
end

% Initialize the matrix and the vector
A = zeros(NIN,NIN);
b = zeros(NIN,1);

% Loop over all the triangles
for i = 1:length(t)
    tt = t(i,1:3);
    x = p(tt, 1:2);
    
    %     Area = area_triangle(x);
    
    % Next step is to calculate the oc-efficients for the basis functions
    M = [ones(3,1),x]; % van der Monde matrix
    v1 = [1,0,0]'; % right hand side for phi1
    v2 = [0,1,0]'; % right hand side for phi2
    v3 = [0,0,1]'; % right hand side for phi3
    % Solve
    phi = M\[v1 v2 v3];
    
    % Gradient of the basis functions
    gphi = phi(2:3,:)';
    
    Area = 0.5*det(M);
    
    for j = 1:3
        for k = 1:3
            if(tt(j) <= NIN && tt(k) <= NIN)
                A(tt(j), tt(k)) = A(tt(j), tt(k)) + Area * gphi(j, :) * (gphi(k, :))';
            end
        end
        if(tt(j) <= NIN)
            b(tt(j)) = b(tt(j)) + (1/3) * Area * f(p(tt(j),1:2));
        end
    end
end

% u = pcg(A,b);
u=A\b;
u = [u; zeros(length(p)-NIN,1)];
plotsolution(p,t,u);

function [g] = f(p)
g = 1;
end

function Area = area_triangle(x)
a1 = edgelength(x([1,2],:));
a2 = edgelength(x([2,3],:));
a3 = edgelength(x([3,1],:));
s = (a1+a2+a3)/2;

Area = sqrt(s*(s-a1)*(s-a2)*(s-a3));
end

function l = edgelength(X)
p1 = X(1,:);
p2 = X(2,:);
dp = p1-p2;
l = sqrt(dp(1).^2 + dp(2).^2);
end

function plotsolution(p,t,u)
figure;
colormap jet;
trisurf(t,p(:,1),p(:,2),u,'EdgeColor','none');
shading interp;
view (0,90);
colorbar;
grid off;
axis tight;
end
