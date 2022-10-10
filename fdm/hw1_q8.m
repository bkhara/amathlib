m = 8;

x0 = 0;
xl = 10;

h = (xl-x0)/(m-1);
A = zeros(m,m);
b = zeros(m,1);

af = [-25/12 4 -3 4/3 -1/4]./h;
A(1,1:5) = af;
A(1,1) = A(1,1)-1;

ab = [1/4 -4/3 3 -4 25/12]./h;
A(m,m-4:m) = ab;
A(m,m) = A(m,m) +1;

cf = [5/6 -5/4 -1/3 7/16 -1/2 1/12]./h^2;
cb = [1/12 -1/2 7/16 -1/3 5/4 -5/6]./h^2;

A(2,1:6) = cf;
A(m,m-5:m) = cb;

%%%%%%%% f-vector
x = linspace(x0,xl,m);
f = -exp(x').*h^2;

d = [-1/12 4/3 -5/2 4/3 -1/12];
for i = 3:(m-2)
    A(i,i-2:i+2) = d;
end
Aleft = zeros(3,6);
Aright = zeros(3,6);

Aleft(1,1:5) = w1;
Aleft(1,3) = -1;
Aleft(2,1:5) = w2;
Aleft(3,2:6) = w2;

Aright(1,2:6) = w1;
Aright(1,4) = 1;
Aright(2,2:6) = w2;
Aright(3,1:5) = w2;

bleft = [0 f(1) f(2)]'.*h^2;
bright = [0 f(m-1) f(m)]'.*h^2;

Bl = [Aleft,bleft];
Br = [Aright,bright];

Bl(2,:) = Bl(2,:) - Bl(1,:).*(Bl(2,1)/Bl(1,1));
Bl(3,:) = Bl(3,:) - Bl(2,:).*(Bl(3,2)/Bl(2,2));

Br(2,:) = Br(2,:) - Br(1,:).*(Br(2,6)/Br(1,6));
Br(3,:) = Br(3,:) - Br(2,:).*(Br(3,5)/Br(2,5));

A(1,1:4) = Bl(3,3:6);
A(m,m-3:m) = Br(3,1:4);

b(1) = Bl(3,end);
b(3:m-2) = f(3:m-2).*h^2;
b(m) = Br(3,end);
%########### b-vector ###################
u = A\b;

uan = exp(10)/2/cos(10).*sin(x) - 0.5.*exp(x);

norm(u-uan)