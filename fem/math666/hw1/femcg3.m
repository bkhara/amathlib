clear;
close all;

% Check with the standard element matrix for beam
% h=1;
% Ke = [12/h, 6, -12/h, 6; 6, 4*h, -6, 2*h; -12/h, -6, 12/h, -6; 6, 2*h, -6, 4*h]./h.^2;
% Ke = [12 -12 6 6;-12 12 -6 -6;6 -6 4 2;6 -6 2 4];

nelV = [10 20 50 100 150 200];
numcases = length(nelV);
hV = 1./nelV;
l2errV = zeros(1,numcases);
enormV = zeros(1,numcases);

X0 = 0;
X1 = 1;
for icases = 1:numcases
    nel = nelV(icases);
    disp(nel);
    nnodes = nel+1;
    
    Xg = linspace(X0, X1, nel+1); % linspace will create nel subdivisions
    
    Ag = zeros(2*nnodes, 2*nnodes);
    bg = zeros(2*nnodes, 1);
    
    % Element assembly starts
    for iel = 1:nel
        xv = Xg([iel, iel+1]);        
        [Ae, be] = ElementMatrixVectorConstruction(xv);
        [Ag, bg] = AssembleElement(Ag, bg, Ae, be, iel);
    end
    
    % Apply boundary condition
    [Ag1, bg1] = ApplyBoundaryConditions(Ag, bg);

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
nu = size(Xg,2);
plot(Xg,u(1:nu),'r', 'LineWidth', 2);
pi2 = pi^2;
uan = @(x) (sin(2*pi.*x).*((18+2.*pi2).*x.^2 + (-24-4*pi2).*x.^3 + (9+2*pi2).*x.^4));
% uan = @(X)((X.^4)/24 - (X.^3)/6 + (X.^2)/4);
fplot(uan, [0,1], '*b', 'LineWidth', 2); 
xlabel('x');
ylabel('u, u_h');
legend('u_h (solution)', 'u (analytical)');

logh=-log10(hV);
loge=log10(l2errV);
logen=log10(enormV);
subplot(1,2,2);
plot(logh, loge, 'o-r', logh, logen, '*-b', 'LineWidth', 2);
xlabel('-log(h)');
ylabel('log(e_{L2}), log(e_{EN})');
legend("||(u-u_h)||_{L2}","||(u-u_h)||_{E}");
[pl2,s] = polyfit(logh, loge, 1);
[pen,s] = polyfit(logh, logen, 1);

end

function [l2err, enorm] = CalcErrorNorms(u, Xg, nel)
nunknowns = 2*size(Xg,2);
l2err = 0;
enorm = 0;
for iel = 1:nel
    i1 = iel;
    i2 = nunknowns/2 + iel;

    xv = Xg([iel, iel+1]);     
    x1 = xv(1);
    x2 = xv(2);
    
    uv = u([i1, i1+1]);
    u1 = uv(1);
    u2 = uv(2);

    duv = u([i2, i2+1]);
    du1 = duv(1);
    du2 = duv(2);

    h = x2 - x1;

    N1 = @(x) (((x - 1).^2.*(x + 2))./4);
    N2 = @(x) (-((x + 1).^2.*(x - 2))./4);
    N3 = @(x) ((h.*(x - 1).^2.*(x + 1))./8);
    N4 = @(x) ((h.*(x - 1).*(x + 1).^2)./8);
    dN1 = @(x) ((3.*x)./2);
    dN2 = @(x) (-(3.*x)./2);
    dN3 = @(x) ((h.*(2.*x - 2))./4 + (h.*(x + 1))./4);
    dN4 = @(x) ((h.*(2.*x + 2))./4 + (h.*(x - 1))./4);

    pi2 = pi^2;
    uan = @(x) (sin(2*pi.*x).*((18+2.*pi2).*x.^2 + (-24-4*pi2).*x.^3 + (9+2*pi2).*x.^4));
    duan = @(x) (sin(2.*pi.*x).*(4.*pi2 - 6.*x.*(4.*pi2 + 24) + 12.*x.^2.*(2.*pi2 + 9) + 36) - 4.*pi.^2.*sin(2.*pi.*x).*((2.*pi2 + 9).*x.^4 + (- 4.*pi2 - 24).*x.^3 + (2.*pi2 + 18).*x.^2) + 4.*pi.*cos(2.*pi.*x).*(2.*x.*(2.*pi2 + 18) + 4.*x.^3.*(2.*pi2 + 9) - 3.*x.^2.*(4.*pi2 + 24)));
    
    usol = @(x)(N1(x).*u1+N2(x).*u2+N3(x).*du1+N4(x).*du2);
    dusol = @(x)(dN1(x).*u1+dN2(x).*u2+dN3(x).*du1+dN4(x).*du2).*(2/h).^2;
    
    l2errel = integral(@(x)(usol(x) - uan(((x1+x2) + h.*x)./2)).^2, -1, 1);
    enormel = integral(@(x)(dusol(x) - duan(((x1+x2) + h.*x)./2)).^2, -1, 1);
    
    J = h/2;
    
    l2err = l2err+l2errel*J;
    enorm = enorm+enormel*J;
end
l2err = sqrt(l2err);
enorm = sqrt(enorm);
end
function [Ae, be] = ElementMatrixVectorConstruction(xv)
Ae = zeros(4,4);
be = zeros(4,1);

x1 = xv(1);
x2 = xv(2);
h = x2 - x1;

N1 = @(x) (((x - 1).^2.*(x + 2))./4);
N2 = @(x) (-((x + 1).^2.*(x - 2))./4);
N3 = @(x) ((h.*(x - 1).^2.*(x + 1))./8);
N4 = @(x) ((h.*(x - 1).*(x + 1).^2)./8);
dN1 = @(x) ((3.*x)./2);
dN2 = @(x) (-(3.*x)./2);
dN3 = @(x) ((h.*(2.*x - 2))./4 + (h.*(x + 1))./4);
dN4 = @(x) ((h.*(2.*x + 2))./4 + (h.*(x - 1))./4);

pi2 = pi^2;
pi4 = pi^4;
pi6 = pi^6;

p1 = @(x) (-64.*pi.*(18 - 6*pi4.*x.^2 + 18.*pi2.*x.^3 + 2.*pi4.*x - 36.*pi2.*x.^2 -27.*x + 12.*pi2.*x + 4.*pi4.*x.^3 + 3.*pi2));
p2 = @(x) (216 - 816.*pi2 - 96.*pi4 + 32.*pi6.*x.^2 - 384.*pi4.*x.^3  -64.*pi6.*x.^3 -288.*pi4.*x.^2 + 3456.*pi2.*x + 576.*pi4.*x - 2592.*pi2.*x.^2 + 144.*pi4.*x.^4 + 32.*pi6.*x.^4);
f = @(x) (p1(x).*cos(2.*pi.*x) + p2(x).*sin(2.*pi.*x));

Ae(1,1) = integral(@(x)(dN1(x).*dN1(x).*(2/h).^4), -1, 1);
Ae(1,2) = integral(@(x)(dN1(x).*dN2(x).*(2/h).^4), -1, 1);
Ae(1,3) = integral(@(x)(dN1(x).*dN3(x).*(2/h).^4), -1, 1);
Ae(1,4) = integral(@(x)(dN1(x).*dN4(x).*(2/h).^4), -1, 1);

Ae(2,2) = integral(@(x)(dN2(x).*dN2(x).*(2/h).^4), -1, 1);
Ae(2,3) = integral(@(x)(dN2(x).*dN3(x).*(2/h).^4), -1, 1);
Ae(2,4) = integral(@(x)(dN2(x).*dN4(x).*(2/h).^4), -1, 1);

Ae(3,3) = integral(@(x)(dN3(x).*dN3(x).*(2/h).^4), -1, 1);
Ae(3,4) = integral(@(x)(dN3(x).*dN4(x).*(2/h).^4), -1, 1);

Ae(4,4) = integral(@(x)(dN4(x).*dN4(x).*(2/h).^4), -1, 1);

Ae(2,1) = Ae(1,2);
Ae(3,1) = Ae(1,3);
Ae(3,2) = Ae(2,3);
Ae(4,1) = Ae(1,4);
Ae(4,2) = Ae(2,4);
Ae(4,3) = Ae(3,4);



be(1) = integral(@(x)(N1(x).*f(((x1+x2) + h.*x)./2)), -1, 1);
be(2) = integral(@(x)(N2(x).*f(((x1+x2) + h.*x)./2)), -1, 1);
be(3) = integral(@(x)(N3(x).*f(((x1+x2) + h.*x)./2)), -1, 1);
be(4) = integral(@(x)(N4(x).*f(((x1+x2) + h.*x)./2)), -1, 1);

J = h/2;

coeff = J;

Ae = Ae.*coeff;
be = be.*coeff;
end

function [Ag, bg] = AssembleElement(Ag, bg, Ae, be, iel)
nunknowns = size(Ag,1);

i1 = iel;
i2 = nunknowns/2 + iel;

Ag(i1:i1+1,i1:i1+1) = Ag(i1:i1+1,i1:i1+1) + Ae(1:2,1:2);
Ag(i1:i1+1,i2:i2+1) = Ag(i1:i1+1,i2:i2+1) + Ae(1:2,3:4);

Ag(i2:i2+1,i1:i1+1) = Ag(i2:i2+1,i1:i1+1) + Ae(3:4,1:2);
Ag(i2:i2+1,i2:i2+1) = Ag(i2:i2+1,i2:i2+1) + Ae(3:4,3:4);

bg(i1:i1+1,1) = bg(i1:i1+1,1) + be(1:2,1);
bg(i2:i2+1,1) = bg(i2:i2+1,1) + be(3:4,1);

end

function [Ag, bg] = ApplyBoundaryConditions(Ag, bg)
nunknowns = size(Ag,1);
Ag(1,:) = 0;
Ag(nunknowns/2 + 1, :) = 0;

Ag(1,1) = 1;
Ag(nunknowns/2 + 1, nunknowns/2 + 1) = 1;

bg(1,1) = 0;
bg(nunknowns/2 + 1, 1) = 0;

end