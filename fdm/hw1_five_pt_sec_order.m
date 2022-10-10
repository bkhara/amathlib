syms h;
assume(h,'real');
A = sym('A',[5,5]);

b = sym('b',[5,1]);

b(1)=0;
b(2)=0;
b(3)=2/h.^2;
b(4)=0;
b(5)=0;
