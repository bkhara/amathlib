clear
close all
nelV = [10 20 50 100 200 500];%, 32, 100, 1000];
numcases = length(nelV);
hV = 2.*pi./nelV;
l2errV = zeros(1,numcases);
enormV = zeros(1,numcases);
X0 = -pi;
X1 = pi;
for icases = 1:numcases
    nel = nelV(icases);
    disp(nel);
    nnodes = nel+1;
    
    Xg = linspace(X0, X1, nel+1); % linspace will create nel subdivisions
    
    Ag = zeros(nnodes, nnodes);
    bg = zeros(nnodes, 1);
    
    % Element assembly starts
    for iel = 1:nel
        xv = Xg([iel, iel+1]);        
        [Ae, be] = ElementMatrixVectorConstruction(xv);
        [Ag, bg] = AssembleElement(Ag, bg, Ae, be, iel);
    end
    
    % Apply boundary condition
    [Ag1, bg1] = ApplyPeriodicBoundaryConditions(Ag, bg);

    % Solve
    u = Ag1\bg1;
    % Calculate error norms
    [l2errV(icases),enormV(icases)] = CalcErrorNorms(u,Xg,nel);
end
% solution plot
% convergence plot
[pl2, pen] = PlotSolutionAndErrorNorms(Xg, u, hV, l2errV, enormV);

function [pl2, pen] = PlotSolutionAndErrorNorms(Xg, u, hV, l2errV, enormV)
% figure;
fig = figure('DefaultAxesFontSize',12,'units','normalized','outerposition',[0 0 0.45 0.35]);
% Plot solution
subplot(1,2,1);
hold on;
plot(Xg,u,'r', 'LineWidth', 2);
uan = @(x)(exp(sin(x)).*exp(cos(x)));
fplot(uan, [-pi,pi], '*b', 'LineWidth', 2); 
xlabel('x');
ylabel('u, u_h');
legend('u_h (solution)', 'u (analytical)','Location','northwest');

% Calculate norms
logh=-log10(hV);
loge=log10(l2errV);
logen = log10(enormV);
[pl2,s] = polyfit(logh, loge, 1);
[pen,s] = polyfit(logh, logen, 1);

% plot error norms
subplot(1,2,2);
plot(logh, loge, 'o-r', logh, logen, '*-b', 'LineWidth', 2);
xlabel('-log(h)');
ylabel('log(e_{L2}), log(e_{EN})');
legend("||(u-u_h)||_{L2}","||(u-u_h)||_{E}");
end

function [l2err, enorm] = CalcErrorNorms(u, Xg, nel)
l2err = 0;
enorm = 0;
for iel = 1:nel
    xv = Xg([iel, iel+1]);     
    x1 = xv(1);
    x2 = xv(2);
    
    uv = u([iel, iel+1]);
    u1 = uv(1);
    u2 = uv(2);

    h = x2 - x1;

    N1 = @(x)((x2-x)/h);
    N2 = @(x)((x-x1)/h);
    dN = [-1,1]./h;

    uan = @(x)(exp(sin(x)).*exp(cos(x)));
    duan = @(x)((cos(x)-sin(x)).*exp(sin(x)).*exp(cos(x)));
    usol = @(x)(N1(x).*u1+N2(x).*u2);
    dusol = @(x)(dN(1).*u1+dN(2).*u2);

    l2errel = integral(@(x)(usol(x) - uan(x)).^2, x1,x2);
    enormel = integral(@(x)(dusol(x) - duan(x)).^2, x1,x2);
    
    l2err = l2err+l2errel;
    enorm = enorm+enormel;
end
l2err = sqrt(l2err);
enorm = sqrt(enorm);
end
function [Ae, be] = ElementMatrixVectorConstruction(xv)
Ae = zeros(2,2);
be = zeros(2,1);

x1 = xv(1);
x2 = xv(2);

h = x2 - x1;

N1 = @(x)((x2-x)/h);
N2 = @(x)((x-x1)/h);
dN = [-1,1]./h;

q = @(x) (3 - sin(x) - sin (2.*x) - cos(x));
f = @(x) (2.*exp(sin(x)).*exp(cos(x)));

Ae(1,1) = integral(@(x)(dN(1)*dN(1) + N1(x).*N1(x).*q(x)), x1, x2);
Ae(1,2) = integral(@(x)(dN(1)*dN(2) + N1(x).*N2(x).*q(x)), x1, x2);
Ae(2,1) = integral(@(x)(dN(2)*dN(1) + N2(x).*N1(x).*q(x)), x1, x2);
Ae(2,2) = integral(@(x)(dN(2)*dN(2) + N2(x).*N2(x).*q(x)), x1, x2);

be(1) = integral(@(x)(N1(x).*f(x)), x1, x2);
be(2) = integral(@(x)(N2(x).*f(x)), x1, x2);
end

function [Ag, bg] = AssembleElement(Ag, bg, Ae, be, iel)
nbf = size(Ae,1);
Ag(iel:iel+nbf-1, iel:iel+nbf-1) = Ag(iel:iel+nbf-1, iel:iel+nbf-1) + Ae;
bg(iel:iel+nbf-1, 1) = bg(iel:iel+nbf-1, 1) + be;
end

function [Ag, bg] = ApplyPeriodicBoundaryConditions(Ag, bg)
nnodes = size(Ag,1);
Ag(nnodes,:) = Ag(nnodes,:) + Ag(1,:); % add the first row to the last row (both the periodic boundaries)
bg(nnodes,1) = bg(nnodes,1) + bg(1,1); % same for bg

Ag(1,:) = 0; % make the first row zero
Ag(1,1) = 1; % introduce the boundary condition
Ag(1,nnodes) = -1; % introduce the boundary condition
bg(1,1) = 0; % basically u_1 = u_n = 0
end