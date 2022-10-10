syms x h N1 N2 M1 M2 uan duan pi2;

N1 = 0.25.*(1-x).^2.*(2+x)
N2 = 0.25.*(1+x).^2.*(2-x)

M1 = h.*(1-x).^2.*(x+1)./8
M2 = h.*(1+x).^2.*(x-1)./8

dN1 = diff(N1,x,2)
dN2 = diff(N2,x,2)

dM1 = diff(M1,x,1)
dM2 = diff(M2,x,1)

uan = (sin(2*pi.*x).*((18+2.*pi2).*x.^2 + (-24-4*pi2).*x.^3 + (9+2*pi2).*x.^4));
duan = diff(uan, x, 2)

% disp(N1)
% disp(N2)
% disp(M1)
% disp(M2)
% 
% disp(dN1)
% disp(dN2)
% disp(dM1)
% disp(dM2)
