%   Example of mesh generation -- square with hole mesh
function [p,t,NIN]=sample_mesh(h)

disp([' ']);
disp(['Square with hole, h=',num2str(h)]);
fdx=inline('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.4))','p');
fh=inline('min(2*sqrt(sum(p.^2,2))+1/2,2)','p');
box=[-1,-1;1,1];
fix=[-1,-1;-1,1;1,-1;1,1];
%[p,t]=distmesh2d(fdx,@huniform,h,box,fix);
[p,t]=distmesh2d(fdx,fh,h,box,fix);
[p,t]=fixmesh(p,t);
post(p,t,@huniform)

[p,t,NIN]=boundary_reorder(p,t,h,fdx);

figure(1); hold on; plot(p(NIN+1:end,1),p(NIN+1:end,2),'ro'); hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post(p,t,fh,varargin)

q=simpqual(p,t);
u=uniformity(p,t,fh,varargin{:});
disp(sprintf(' - Min quality %.2f',min(q)))
disp(sprintf(' - Uniformity %.1f%%',100*u))
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%