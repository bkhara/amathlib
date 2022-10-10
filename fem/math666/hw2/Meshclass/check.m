
clear all;
clc;
close all;

% fouter=@(p) sqrt(p(:,1).^2+p(:,2).^2)-0.75-0.25*sin(5*atan(p(:,2)./p(:,1)));
% finner= @(p) 0.25*0.25-p(:,1).^2-p(:,2).^2;
% % pfix=[-1,-1;-1,1;1,-1;1,1]
% %  [p,t]=distmesh2d(fouter,@huniform,0.01,[-1,-1;1,1],[]);
% 
% % fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-0.5;
% %     [p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
% 

counter = 1;
for theta =0:0.02:2*pi
    r(counter) = 0.75+0.25*sin(5*theta);
    Th(counter) = theta;
    counter = counter + 1;
end
% polarplot(Th,r);
X = r.*cos(Th);
Y = r.*sin(Th);
pfix=[X(1:end); Y(1:end)]';

Th_circ = 0:0.1:2*pi;
X = 0.25*cos(Th_circ);
Y = 0.25*sin(Th_circ);
pfix2 = [X;Y]';
pfix = vertcat(pfix,pfix2);

% fdx=inline('ddiff(dstar(p,1),dcircle(p,0,0,0.5))','p');
fdx=inline('ddiff(dstar(p),dcircle(p,0,0,0.25))','p');
% fd = @(p) dstar(p);
[p,t]=distmesh2d(fdx,@huniform,0.05,[-1,-1;1,1],pfix);

