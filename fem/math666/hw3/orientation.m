function [tri_orient, new_order] = orientation(X,Y,nodelist)
xp = X(nodelist);
yp = Y(nodelist);
new_order = zeros(1,3);
tri_orient = 1; % upright

ymin = min(yp);
ymax = max(yp);

lur = ymax == yp;
ldr = ymin == yp;

if(sum(lur) == 1)
    tri_orient = 1;
elseif(sum(ldr) == 1)
    tri_orient = 0; % downright
end

if (tri_orient) % upright
    new_order(3) = nodelist(yp==ymax);
    temp = xp(yp~=ymax);
    
    lmin = xp==min(temp); % will be length 2, has to be resolved
    nd = nodelist(lmin);
    temp2 = nd ~= new_order(3);
    new_order(1) = nd(temp2);
    
    lmax = xp==max(temp); % will be length 1, so no problem
    new_order(2) = nodelist(lmax);    
else
    new_order(1) = nodelist(yp==ymin);
    temp = xp(yp~=ymin);
    
    lmin = xp==min(temp); % will be length 1, so no problem
    new_order(2) = nodelist(lmin); 
    
    lmax = xp==max(temp); % will be length 2, has to be resolved
    nd = nodelist(lmax);
    temp2 = nd ~= new_order(1);
    new_order(3) = nd(temp2);
end
disp(nodelist);
disp(tri_orient);
end