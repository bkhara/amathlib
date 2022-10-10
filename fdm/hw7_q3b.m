% Fourier spectral method for solving the
% 1D variable coefficient wave equation on [0,2*pi]
% with periodic boundary conditions

% Grid, variable coefficient, and initial data:
clear;
close all;
nv = 2.^(5:1:8);
ncases = length(nv);
hv = zeros(ncases,1);
ev = zeros(ncases,1);
for icase = 1:ncases
    N = nv(icase);
    h = 2*pi/N;
    x = h*(1:N);
    t = 0;
    dt = h/4;
    c = .2 + sin(x-1).^2;
    v = exp(-100*(x-1).^2);
    vold = exp(-100*(x-.2*dt-1).^2);
    u0 = v;
    
    tmax = 10*pi/sqrt(6);
    
    nplots = round(tmax/dt);
    for i = 1:nplots
        if(abs(t-tmax)<dt)
            t = tmax;
        else
            t = t+dt;
        end
        v_hat = fft(v);
        w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
        w = real(ifft(w_hat));
        vnew = vold - 2*dt*c.*w;
        vold = v;
        v = vnew;
    end
    vdiff = abs(u0 - v);
    hv(icase) = h;
    ev(icase) = norm(vdiff,inf);
end

[p,s] = polyfit(log10(hv),log10(ev),1);
disp(['Rate of Convergence = ',num2str(p(1))]);

% solution and convergence plto
figure;
hold on;
plot(x,u0,'ro','LineWidth',2);
plot(x,v,'b-','LineWidth',1.5);
legend('u_0','u_T');
bb=xlabel('$u_0,u_T$','FontSize',16);
set(bb,'Interpreter','latex');
box on;

% subplot(1,3,3);
figure;
loglog(hv,ev,'bo-','LineWidth',2);
bb=xlabel('$\log (h)$','FontSize',16);
set(bb,'Interpreter','latex');
bb=ylabel('$\log_{10}\left(\|e\|_{L^{\infty}}\right)$','FontSize',16);
set(bb,'Interpreter','latex');
box on;