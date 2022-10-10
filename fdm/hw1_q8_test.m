clear;
close all
mV = [10 20 50 100 200 500 1000];
eV = zeros(1,length(mV));
hV = zeros(1,length(mV));

for im = 1:length(mV)
    m = mV(im);
    
    x0 = 0;
    xl = 10;
    
    h = (xl-x0)/(m-1);
    hV(im) = h;
    
    A = zeros(m,m);
    b = zeros(m,1);
    
    %%%%%%%% A-matrix
    % first row
    af = [-25 48 -36 16 -3]./12;
    A(1,1:5) = af./h;    
    % last row
    ab = [3 -16 36 -48 25]./12;
    A(m,m-4:m) = ab./h;    
    % Boundary conditions
    A(1,1) = A(1,1)-1; % left end
    A(m,m) = A(m,m) +1; % right end
    
    % the second row and the second-last row
    cf = [10 -15 -4 14 -6 1]./12;
    cb = [1 -6 14 -4 -15 10]./12;    
    A(2,1:6) = cf./h^2;
    A(m-1,m-5:m) = cb./h^2;
    
    % all the other interior rows
    d = [-1 16 -30 16 -1]./12;
    d = d./h^2;
    for i = 3:(m-2)
        A(i,i-2:i+2) = d;
    end
    
    % Adding 1 to each diagonal of the interior node to account 
    % for the linear u-term in the BVP: u'' + u = f
    for i = 2:m-1
        A(i,i) = A(i,i) + 1;
    end
    
    %%%%%%%% f-vector
    x = linspace(x0,xl,m);
    f = -exp(x');
    f(1) = 0;
    f(m) = 0;
    
    u=A\f;
    uan1 = (exp(10)/2/cos(10)).*(sin(x) + cos(x)) - 0.5.*exp(x);
    uan1 = uan1';
    
    eV(im) = norm(f-A*uan1);
end
[p,s] = polyfit(log10(hV),log10(eV),1)
figure;
hold on;
plot(x,u,'r-');
x = linspace(x0,xl,floor(mV(end)/4));
uan = (exp(10)/2/cos(10)).*(sin(x) + cos(x)) - 0.5.*exp(x);
plot(x,uan','b*');
legend('u','uan');
xlabel('x');
ylabel('u(x)');

figure;
plot(log10(hV),log10(eV),'bo-')
aa=xlabel('$\log_{10}(h)$','FontSize',16);
set(aa,'Interpreter','latex');
bb=ylabel('$\log_{10}\left(\|f-Au_{exact}\|\right)$','FontSize',16);
set(bb,'Interpreter','latex');
title('Convergence of the FDM solution of the BVP');