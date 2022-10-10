close all;
ndiv = 100;

h = 1e-1;
klim = 20;
klimits = [-klim,klim];
k = (linspace(klimits(1),klimits(2),ndiv))';

ginf = @(k,h) (k);
g2 = @(k,h) (sin(k.*h)./h);
g4 = @(k,h) ((8.*sin(k.*h)-sin(2.*k.*h))./(6*h));

ginfp = @(k) (ginf(k,h));
g2p = @(k) (g2(k,h));
g4p = @(k) (g4(k,h));

figure;
hold on;

plot(k,ginfp(k),'r','LineWidth',2);
plot(k,g2p(k),'b','LineWidth',2);
plot(k,g4p(k),'k','LineWidth',2);
legend('g_{\infty}','g_{2}','g_{4}','Location','northwest');
xlabel('k');
box on;