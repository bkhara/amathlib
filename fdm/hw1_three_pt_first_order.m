syms h;
assume(h,'real');

npt = 3;

A = sym('A',[npt,npt]);

A(1,1)=1;
A(1,2)=1;
A(1,3)=1;

A(2,1)=0;
A(2,2)=h;
A(2,3)=2.*h;

A(3,1)=0;
A(3,2)=0.5.*h.^2;
A(3,3)=2.*h.^2;

b1 = sym('b1',[npt,1]);
b1(1)=0;
b1(2)=1;
b1(3)=0;

b2 = sym('b2',[npt,1]);
b2(1)=0;
b2(2)=-1;
b2(3)=0;
