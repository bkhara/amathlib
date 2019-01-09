clc;
clear;
m = 5;
imax=11;
a = randi(imax,m);
% a=magic(m);
a = 0.5.*(a+a');

% ############################################################
% TESTING THE hessenberg AND THE schurQRbasic ROUTINES

% [p,h] = hessenberg(a); % defined as A=PHP'
% [q,t] = schurQRbasic(h); % deefined as H=QTQ'
% u = p*q; % u is defined as A=UTU'=(PQ)T(PQ)'
% 
% [va,da]=eig(a);
% [vh,dh]=eig(h);
% 
% das=sort(diag(da));
% dhs=sort(diag(dh));
% dts=sort(diag(t));
% % disp(sort(da)-sort(dh));
% 
% schur = q*t*q';
% 
% % dts-das
% norm(dts-das,1)
% norm(dts-das,2)

% ############################################################
% TESTING THE eigen ROUTINE

[va,da]=eig(a);
[vm,dm]=eigen(a);