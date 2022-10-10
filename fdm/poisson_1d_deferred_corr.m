% This code works for the 1-D BVP
% -u_xx = f
% u = g on boundary (fully Dirichlet conditions)
clear;
close all
mV = [10 20 50 100 200 500];
eV = zeros(1,length(mV));
hV = zeros(1,length(mV));

for im = 1:length(mV)
    m = mV(im);
    
    x0 = 0;
    xl = 1;
    
    h = (xl-x0)/(m-1);
    hV(im) = h;
    
    A = zeros(m,m);
    b = zeros(m,1);
    
    %%%%%%%% A-matrix
    % first row
    A(1,1) = 1;
    % last row
    A(m,m) = 1;
    % all the other interior rows
    d = [-1 2 -1]./h^2;
    for i = 2:(m-1)
        A(i,i-1:i+1) = d;
    end
    
    %%%%%%%% f-vector
    x = (linspace(x0,xl,m))';
    f = -fEval(x);
    f(1) = uBC(x0);
    f(m) = uBC(xl);
    
    u=A\f;
    uan1 = uAnalytic(x);%(exp(10)/2/cos(10)).*(sin(x) + cos(x)) - 0.5.*exp(x);
    %uan1 = uan1';
    
    eV(im) = norm(f-A*uan1,inf);
end
[p,s] = polyfit(log10(hV),log10(eV),1)
figure;
hold on;
plot(x,u,'r-');
x = (linspace(x0,xl,floor(mV(end)/4)))';
uan = uAnalytic(x);%(exp(10)/2/cos(10)).*(sin(x) + cos(x)) - 0.5.*exp(x);
plot(x,uan,'b*');
legend('u','uan');
xlabel('x');
ylabel('u(x)');

figure;
plot(log10(hV),log10(eV),'bo-')
aa=xlabel('$\log_{10}(h)$','FontSize',16);
set(aa,'Interpreter','latex');
bb=ylabel('$\log_{10}\left(\|f-Au_{exact}\|\right)$','FontSize',16);
set(bb,'Interpreter','latex');
title('Convergence of the FDM solution using 2nd order Central Difference');


function [u] = uAnalytic(x)
u = exp(x);
end
function [f] = fEval(x)
f = exp(x);
end
function [u] = uBC(x)
u = uAnalytic(x);
end