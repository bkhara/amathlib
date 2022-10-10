% This code works for the 1-D BVP
% -u" = f
% u(0) = alpha, u'(L) = sigma
clear;
close all;
global icase;
icase = 5;
mV = [10 20 50 100 200 500 1000];
eV = zeros(1,length(mV));
hV = zeros(1,length(mV));

for im = 1:length(mV)
    m = mV(im);
    
    x0 = 0;
    xl = 1;
    
    h = (xl-x0)/(m-1);
    hV(im) = h;
    
    x = linspace(x0,xl,m);
    x = x';
    A = zeros(m,m);
    b = zeros(m,1);
    f = zeros(m,1);
    uan = zeros(m,1);
    
    bc1 = EvalFunc(x0);
    bc2 = (h/2)*EvalForce(xl) + EvalDer(xl);
    
    %%%%%%%% A-matrix
    % first row
    A(1,1) = 1;
    % last row
    ab = [-1 1];
    A(m,m-1:m) = ab./h;
    
    % all the other interior rows
    d = [-1 2 -1];
    d = d./h^2;
    for i = 2:(m-1)
        A(i,i-1:i+1) = d;
    end
    
    %%%%%%%% f-vector
    for ii = 1:m
        uan(ii) = EvalFunc(x(ii));
        f(ii) = EvalForce(x(ii));%     eV(im) = norm(F-B\uan,2);
    end
    f(1) = bc1;
    f(m) = bc2;
    
    %%%%%%%% Calculating B=inv(A) using Green's function
    B = zeros(m,m); % B=inv(A)
    for i = 1:m % loop on rows
        B(i,1) = 1;
        B(i,m) = x(i);
    end
    for j = 2:m-1 % loop on columns
        for i = 1:m
            if(i<=j)
                B(i,j) = h*x(i);
            else
                B(i,j) = h*x(j);
            end
        end
    end
    
    uFD = A\f; % solution by linear solve
    uG  = B*f; % solution using the Green's function
    
    % error in L2
    eV(im) = norm(f-A*uan,2);
end
fig=1;
if(fig)
    figure;
    hold on;
    plot(x,uFD,'r-');
    plot(x,uG,'ko');
    y = linspace(x0,xl,floor(mV(end)/4));
    uan = EvalFunc(y);
    plot(y,uan','b*');
    legend('u','uan');
    xlabel('x');
    ylabel('u(x)');
end

if(length(mV)>2)
    [p,s] = polyfit(log10(hV),log10(eV),1)
    figure;
    plot(log10(hV),log10(eV),'bo-')
    aa=xlabel('$\log_{10}(h)$','FontSize',16);
    set(aa,'Interpreter','latex');
    bb=ylabel('$\log_{10}\left(\|f-A*u_{exact}\|\right)$','FontSize',16);
    set(bb,'Interpreter','latex');
end

function [u] = EvalFunc(x)
global icase;
switch(icase)
    case 1
        u = exp(x);
    case 2
        u = x.^3;
    case 3
        u = cos(x);
    case 4
        u = sin(x);
    case 5
        u = x.^5;
end
end

function [Du] = EvalDer(x)
global icase;
switch(icase)
    case 1
        Du = exp(x);
    case 2
        Du = 3.*x.^2;
    case 3
        Du = -sin(x);
    case 4
        Du = cos(x);
    case 5
        Du = 5.*x.^4;
end
end

function [f] = EvalForce(x)
global icase;
switch(icase)
    case 1
        f = exp(x);
    case 2
        f = 6.*x;
    case 3
        f = -cos(x);
    case 4
        f = -sin(x);
    case 5
        f = 20.*x.^3;
end
f = -f;
end
