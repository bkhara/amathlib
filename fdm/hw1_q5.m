syms h h1 h2 h3;
assume(h, 'real');
assume(h1, 'real');
assume(h2, 'real');
assume(h3, 'real');

A = [1 1 1 1;...
    -h1     0 h2     (h2+h3);...
     h1^2 0 h2^2 (h2+h3)^2;...
    -h1^3 0 h2^3 (h2+h3)^3];

b = 2.*[0 0 1 0]';
x = A\b;