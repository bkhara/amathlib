n = 1000;
p = linspace(-16,1,n);

hV = 10.^p;
errV = zeros(1,n);

icase = 9;

switch(icase)
    %############ FIRST DERIVATIVES ##############
    case 1
        der = 1; %order=2
        st = [-1 1];
        C = [-1/2 1/2];
    case 2
        der = 1; % order=4
        st = [0 1 2 3 4];
        C = [-25 48 -36 16 -3]./12;
    case 3
        der = 1;
        st = [0 1 2 3 4 5]; % order=5
        C = [-137 300 -300 200 -75 12]./60;
        
        %############ SECOND DERIVATIVES ##############
    case 4
        der = 2; % order=4
        st = [-2 -1 0 1 2];
        C = [-1 16 -30 16 -1]./12;
    case 5
        der = 2; % order=4
        st = [0 1 2 3 4 5];
        C = [45 -154 214 -156 61 -10]./12;
    case 6
        der = 2; % order=4
        st = [-1 0 1 2 3 4];
        C = [10 -15 -4 14 -6 1]./12;
    case 7
        der = 2; % order=4
        st = [-4 -3 -2 -1 0 1];
        C = [1 -6 14 -4 -15 10]./12;
    case 8
        der = 2; % order =
        st = [-4 -3 -2 -1 0];
        C = [11 -56 114 -104 35]./12;
    case 9
        der = 2; % order =
        st = [-1 0 1 2];
        C = [1 -2 1 0];
end

npt = length(st);

x0 = 1;
for i = 1:n
    h = hV(i);
    if(der == 1)
        div = h;
    elseif(der == 2)
        div = h^2;
    end
    coef = C/div;
    
    difsum = 0;
    for j = 1:npt
        difsum = difsum + coef(j) * funcEval(x0 + st(j)*h);
    end
    
    d2_anal = funcEvalD2(x0);
    %     d2_anal = funcEvalD1(x0);
    d2_fd = difsum;
    
    errV(i) = abs(d2_fd - d2_anal);
end
[p,s] = polyfit(log10(hV),log10(errV),1);
slope = p(1);
figure;
hold on;
plot(log10(hV),log10(errV));

if(icase == 4)
    xline1 = linspace(-15,-4,5);
    yline1 = -2.*xline1 -13;
    plot(xline1, yline1, 'r--');
    
    xline2 = linspace(-2,0.5,5);
    yline2 = 4.*xline2-3;
    plot(xline2,yline2,'b--');
    
    xc1 = -2.2;
    xline3 = [1 1 1]'.*xc1;
    yline3 = linspace(-15,20,3);
    plot(xline3,yline3,'k--')
    
    xlabel('log_{10}(h)');
    ylabel('log_{10}(e)');
end




function u = funcEval(x)
u = exp(x);
end
function du = funcEvalD1(x)
du = exp(x);
end
function d2u = funcEvalD2(x)
d2u = exp(x);
end

