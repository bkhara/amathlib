N = 1001;
H = 10.^linspace(-5,-1,N);
errV = zeros(1,N);
err_uniformV = zeros(1,N);

w = zeros(1,4);
x0 = 1;
for i = 1:N
    h1 = H(i)*rand(1);
    h2 = H(i)*rand(1);
    h3 = H(i)*rand(1);
    
    w(1) = (2*(2*h2 + h3))/(h1*(h1 + h2)*(h1 + h2 + h3));
    w(2) =  -(2*(2*h2 - h1 + h3))/(h1*h2*(h2 + h3));
    w(3) = (2*(h2 - h1 + h3))/(h3*(h2^2 + h1*h2));
    w(4) =(2*(h1 - h2))/(h2^2*h3 + 2*h2*h3^2 + h1*h2*h3 + h3^3 + h1*h3^2);
    
    points = [x0-h1, x0, x0+h2, x0+(h2+h3)]';
    uVector = exp(points);
    d2u_fd = w*uVector;
    d2u_anal = exp(1);
    
    errV(i) = abs(d2u_fd-d2u_anal);
end
for i = 1:N
    h1 = H(i);
    h2=h1;h3=h1;
    
    w(1) = (2*(2*h2 + h3))/(h1*(h1 + h2)*(h1 + h2 + h3));
    w(2) =  -(2*(2*h2 - h1 + h3))/(h1*h2*(h2 + h3));
    w(3) = (2*(h2 - h1 + h3))/(h3*(h2^2 + h1*h2));
    w(4) =(2*(h1 - h2))/(h2^2*h3 + 2*h2*h3^2 + h1*h2*h3 + h3^3 + h1*h3^2);
    
    points = [x0-h1, x0, x0+h2, x0+(h2+h3)]';
    uVector = exp(points);
    d2u_fd = w*uVector;
    d2u_anal = exp(1);
    
    err_uniformV(i) = abs(d2u_fd-d2u_anal);
end

lh = log10(H);
le = log10(errV);
leu = log10(err_uniformV);

figure;
hold on;
plot(lh,le);
plot(lh,leu);
legend('Random values of h_1, h_2, h_3 in [0,H(i)]', 'Uniform h_1=h_2=h_3=H(i)');
xlabel('log_{10}(H_i)');
ylabel('log_{10}(e)');

n = 500;
A = [ones(n,1), log10(H(1:n)')];
b = log10(errV(1:n)');
xx = A\b;
p = xx(2)