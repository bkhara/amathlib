syms h h1 h2 h3;
assume(h,'real');
assume(h1,'real');
assume(h2,'real');
assume(h3,'real');


A = [1 1 1 1;...
    -h1 0 h2 (h2+h3);...
    h1.^2 0 h2.^2 (h2+h3).^2;...
    -h1.^3 0 h2.^3 (h2+h3).^3];
b = [0 0 2 0]';

w = A\b;

%     -(2*(2*h2 + h3))/(h1*(h1^2 - 2*h1*h2 - h3*h1 + h2^2 + h3*h2))
%                            (2*(h1 + 2*h2 + h3))/(h1*h2*(h2 + h3))
%                         -(2*(h1 + h2 + h3))/(h3*(- h2^2 + h1*h2))
%  -(2*(h1 + h2))/(h2^2*h3 + 2*h2*h3^2 - h1*h2*h3 + h3^3 - h1*h3^2)


% A = [1 1 1 1;...
%     -h 0 h (h+h);...
%     h.^2 0 h.^2 (h+h).^2;...
%     -h.^3 0 h.^3 (h+h).^3];
% b = [0 0 2 0]';
% 
% w = A\b;