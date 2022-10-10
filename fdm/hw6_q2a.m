syms h x u1 u2 u3 u4;
assume(h,'real');
assume(x,'real');
assume(u1,'real');
assume(u2,'real');
assume(u3,'real');
assume(u4,'real');

A = [1 (x-2*h) (x-2*h)^2 (x-2*h)^3;...
    1 (x-h) (x-h)^2 (x-h)^3;...
    1 x x^2 x^3;...
    1 (x+h) (x+h)^2 (x+h)^3];
    
B = [u1 u2 u3 u4]';

X = A\B;

X = subs(X,x,0)


X=simplify(X)

Ainv = inv(A);
Ainv = subs(Ainv,x,0);
Ainv = simplify(Ainv);

xx = [1 x x^2 x^3];

yy=simplify(xx*Ainv);
yy=yy'
y1=diff(yy,x)
y2=diff(yy,x,2)
y3=diff(yy,x,3)

y10=subs(y1,x,0)
y20=subs(y2,x,0)
y30=subs(y3,x,0)
