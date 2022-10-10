clear;
close all
mV = [10 20 50 100 200 500 1000];
eV = zeros(1,length(mV));
hV = zeros(1,length(mV));

for im = 1:length(mV)
    m = mV(im);
    
    x0 = 0;
    xl = 2*pi;
    
    h = (xl-x0)/(m-1);
    hV(im) = h;
    
    x = linspace(x0,xl,m);
    x = x';
    A = zeros(m,m);
    b = zeros(m,1);
    f = zeros(m,1);
    uan = zeros(m,1);
    
    alpha = EvalFunc(x0);
    sigma = (h/2)*EvalForce(xl) - EvalDer(xl);
    
    %%%%%%%% A-matrix
    % first row
    A(1,1) = 1;
    % last row
    ab = [1 -1];
    A(m,m-1:m) = ab./h;
    
    % all the other interior rows
    d = [1 -2 1];
    d = d./h^2;
    for i = 2:(m-1)
        A(i,i-1:i+1) = d;
    end
    
    %%%%%%%% f-vector
    for ii = 1:m
        uan(ii) = EvalFunc(x(ii));
        f(ii) = EvalForce(x(ii));
    end
    f(1) = alpha;
    f(m) = sigma;
    
    % Calculating inv(A) using Green's function
    B = zeros(m,m); % B=inv(A)
    for i = 1:m % loop on rows
        B(i,1) = 1;
        B(i,m) = -x(i);
    end
    for j = 2:m-1 % loop on columns
        for i = 1:m
            if(i<=j)
                B(i,j) = -h*x(i);
            else
                B(i,j) = -h*x(j);
            end
        end
    end
    
    uFD=A\f;
%     u = B*f;
    
    eV(im) = norm(f-A*uan);
end
fig=1;
if(fig)
    figure;
    hold on;
    plot(x,uFD,'r-');
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
% u = exp(x);
% u = x.^3;
% u = cos(x);
% u = sin(x);
u = x.^5;
end

function [Du] = EvalDer(x)
% Du = exp(x);
% Du = 3.*x.^2;
% Du = -sin(x);
% Du = cos(x);
Du = 5.*x.^4;
end

function [f] = EvalForce(x)
% f = exp(x);
% f = 6.*x;
% f = -cos(x);
% f = -sin(x);
f = 20.*x.^3;
end
