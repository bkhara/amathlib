frac = 0.2489;
h = 1/3;
xp = [0,h,0];
yp = [0,0,h];

% xp = [h,h,0];
% yp = [0,h,h];

x1 = xp(1);
x2 = xp(2);
x3 = xp(3);

y1 = yp(1);
y2 = yp(2);
y3 = yp(3);

M = [1 x1 y1;...
    1 x2 y2;...
    1 x3 y3];

detM = det(M);

Minv = inv(M);

xg = h*frac;
yg = h*frac;

% N1 = sum(Minv(:,1).*[1 xg yg]);
% N2 = sum(Minv(:,2).*[1 xg yg]);
% N3 = sum(Minv(:,3).*[1 xg yg]);

xyg = [1 xg yg];

N1 = dot(Minv(:,1),xyg);
N2 = dot(Minv(:,2),xyg);
N3 = dot(Minv(:,3),xyg);

disp(0.5*detM - 0.5*h.^2)
disp(N1+N2+N3)

E = [(x2*y3 - x3*y2), (x1*y3 - x3*y1), (x1*y2 - x2*y1);...
    (y2 - y3), (y3 - y1), (y1 - y2);...
    (x3 - x2), (x1 - x3), (x2 - x1)];

E = E./det(M);
disp(E)
disp(Minv - E);


