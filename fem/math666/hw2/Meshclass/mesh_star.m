function [p,t,NIN]=mesh_star(h,fh,uniformityflag)
% h = 0.05;
nptstar = 600;
nptcirc = 200;
endflag = 1;

% Star shaped region
tvstar = linspace(0,2*pi,nptstar);
tvstar = tvstar';
rvstar = 0.75 + 0.25 .* sin(5 .* tvstar);
% Boundary fixes for the star shaped region
xstar = rvstar(1:end-1*endflag) .* cos(tvstar(1:end-1*endflag));
ystar = rvstar(1:end-1*endflag) .* sin(tvstar(1:end-1*endflag));
pfixstar = [xstar,ystar];

% Circular hole
tvcirc = linspace(0,2*pi,nptcirc);
tvcirc = tvcirc';
rcirc = 0.25;
% Boundary fixes for the circular hole
xcirc = rcirc .* cos(tvcirc(1:end-1*endflag));
ycirc = rcirc .* sin(tvcirc(1:end-1*endflag));
pfixcirc = [xcirc, ycirc];

% All fixed points
if(uniformityflag == 1)
    pfixall = [pfixstar;pfixcirc];
else
    pfixall = pfixstar;
end


% Distance function handle, that defines the region by the subtraction of
% the outer and the inner boundary
fdx = @(p) (ddiff(dstar(p),dcircle(p,0,0,0.25)));

% Plot the domain
figure;
plot(pfixstar(1:end,1),pfixstar(1:end,2),pfixcirc(1:end,1),pfixcirc(1:end,2));

% Generate the mesh
% fh=@(p) (2 - 2*ddiff(dstar(p),dcircle(p,0,0,0.25)));
if(uniformityflag == 1)
    [p,t]=distmesh2d(fdx,@huniform,h,[-1,-1;1,1],pfixall);
else
    [p,t]=distmesh2d(fdx,fh,h,[-1,-1;1,1],pfixall);
end
size(p)
size(t)
[p,t]=fixmesh(p,t);
size(p)
size(t)
[p,t,NIN]=boundary_reorder(p,t,h,fdx);
figure(1); hold on; plot(p(NIN+1:end,1),p(NIN+1:end,2),'ro'); hold off;
end


