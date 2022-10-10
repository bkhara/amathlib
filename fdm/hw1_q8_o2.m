mV = 10%[10 40 80 160 200 500];
eV = zeros(1,length(mV));

for im = 1:length(mV)
    m = mV(im);
    
    x0 = 0;
    xl = 10;
    
    h = (xl-x0)/(m-1);
    A = zeros(m,m);
    b = zeros(m,1);
    x = linspace(x0,xl,m);
    f = -exp(x');
    
    %     f(1) = 0;
    %     f(m) = 0;
    %     af = [-3 4 -1]./(2*h);
    %     A(1,1:3) = af;
    %     A(1,1) = A(1,1)-1;
    %     ab = [1 -4 3]./(2*h);
    %     A(m,m-2:m) = ab;
    %     A(m,m) = A(m,m)+1;
    
    f(1) = f(1)/2;
    f(m) = f(m)/2;
    A(1,1:2)=[-(1+h) 1]./h^2;
    A(m,m-1:m) = [1 -(1+h)]./h^2;
    
    d = [1 -2 1]./h^2;
    d(2) = d(2)+1;
    for i = 2:(m-1)
        A(i,i-1:i+1) = d;
    end
    
    u=A\f;
    uan1 = (exp(10)/2/cos(10)).*(sin(x) + cos(x)) - 0.5.*exp(x);
    uan1 = uan1';
    %uan2 = -(cos(10)*exp(x) - 2^(1/2)*exp(10)*sin(x + pi/4))/(2*cos(10));
    eV(im) = norm(f-A*uan1)/norm(f);
    
    %     norm(f-A*uan1')/norm(f);
    
end
[p,s] = polyfit(log10(mV),log10(eV),1)

plot(log10(mV),log10(eV))

figure;
hold on;

plot(x,uan1);
plot(x,u);
% plot(x,(f-A*uan1)./max(f));
% legend('u','uan');
legend('uan','u');