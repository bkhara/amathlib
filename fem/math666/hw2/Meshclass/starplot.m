clear;
close all;
h = 0.025;
uniformityflag = 0;
% fh=@(p) (1 - 3 * ddiff(dstar(p),dcircle(p,0,0,0.25)));
fh=@(p) min(2*sqrt(sum(p.^2,2))+0.3,2);
[p,t,NIN]=mesh_star(h, fh, uniformityflag);
post(p,t,@huniform)

save('starmesh','p','t','NIN');


function post(p,t,fh,varargin)

q=simpqual(p,t);
u=uniformity(p,t,fh,varargin{:});
disp(sprintf(' - Min quality %.2f',min(q)))
disp(sprintf(' - Uniformity %.1f%%',100*u))
disp(' ')
end