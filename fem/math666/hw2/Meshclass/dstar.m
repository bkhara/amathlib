function d=dstar(p)
x = p(:,1);
y = p(:,2);
d = sqrt(x.^2 + y.^2) - 0.75 - 0.25 * sin(5*atan2(y,x));
end
