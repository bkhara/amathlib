syms x1 x2 x3;
syms y1 y2 y3;


M = [1 x1 y1;...
    1 x2 y2;...
    1 x3 y3];

pretty(inv(M))