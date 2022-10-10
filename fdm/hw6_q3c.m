ndiv = 30;

x = linspace(-10,1,ndiv);
y = linspace(-50,50,ndiv);
[X,Y]=meshgrid(x,y);
z = abs(exp(X+1i.*Y));
Z = double(z<1);

surf(X,Y,Z)