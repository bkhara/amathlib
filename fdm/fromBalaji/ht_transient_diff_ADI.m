%% ADI in 2D

clear
close all
clc


all_K = [];
all_T_resized = [];
all_T_solution = [];

%%
% thermal diffusitivity
% load('/home/balajip/Desktop/grf_figs/testdata.mat')
% A = h5read('p3_sec_org_0.145_10per_4766.h5', '/images');
load('grf_data.mat');
A = all_morph;
for ind = 1:1
    disp(ind)
    K = reshape(A(:,:, ind), [129,129]);
    K = 9.*imresize(K>0.5, [64,64])+1;
    %     K = ones(size(K));
    [T1, K1, T2] = solve_ht_diffusion(K);
    % save the solutions
    if ind<1000
        S = load('./transient_ht_solutions_1.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_1.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<2000
        S = load('./transient_ht_solutions_2.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_2.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<3000
        S = load('./transient_ht_solutions_3.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_3.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<4000
        S = load('./transient_ht_solutions_4.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_4.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<5000
        S = load('./transient_ht_solutions_5.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_5.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<6000
        S = load('./transient_ht_solutions_6.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_6.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<7000
        S = load('./transient_ht_solutions_7.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_7.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<8000
        S = load('./transient_ht_solutions_8.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_8.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<9000
        S = load('./transient_ht_solutions_9.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_9.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<10000
        S = load('./transient_ht_solutions_10.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_10.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<11000
        S = load('./transient_ht_solutions_11.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_11.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<12000
        S = load('./transient_ht_solutions_12.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_12.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<13000
        S = load('./transient_ht_solutions_13.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_13.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<14000
        S = load('./transient_ht_solutions_14.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_14.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    elseif ind<15000
        S = load('./transient_ht_solutions_15.mat');
        all_K = [S.all_K, {K}];
%         all_T_resized = [S.all_T_resized, {T1}];
        all_T_solution = [S.all_T_solution, {T2}];
        save('./transient_ht_solutions_15.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
        clear S all_K all_T_resized all_T_solution;
    else
        save('./transient_ht_solutions.mat', 'all_K', ...
            'all_T_solution', '-v7.3')
    end
    
end


%% Numerical solution
function [T_resize, K, T_solution] = solve_ht_diffusion(K)
%% Declaration of variables
%%
% Spatial resolution
x_start = 0;
x_end = 1;
y_start = 0;
y_end = 1;
%%
% Number of grid points
N_x = size(K,1);
N_y = size(K,2);
%%
% Element length
del_x = (x_end - x_start)/(N_x - 1);
del_y = (y_end - y_start)/(N_y - 1);
%%
% Positions where the Temperature is calculated
x = linspace(x_start, x_end,N_x);
y = linspace(y_start, y_end,N_y);
%%
alpha = 25;
%%
% Temporal Resolution
t_start = 0;
t_end = 25;
t = t_start;
del_t = alpha*del_x^2;
time = t_start:del_t:t_end;
N_t = length(time);

%%
% Initialisation of Temperature
% T_old represents the temperature calculated at the present time step
T_old = zeros(N_x,N_y);
%%
% Boundary conditions -- one side heated wall
T_old(1,:) = 0;
T_old(N_x,:) = 1;
T_old(:,1) = 0;
T_old(:,N_y) = 0;
iter = 1;

T_solution = zeros([N_x, N_y, N_t]);
T = T_old;
%%
while(t<t_end)
    T_solution(:, :, iter) = reshape(T, [N_x, N_y, 1]);
    % Plotting at specified intervals of time
%     if t == t_start
%         figure
%         hold on
%         contourf(K,20);
%         xlabel('Non- dimensional y-distance');
%         ylabel('Non- dimensional x-distance');
%         str = sprintf('Diffusivity distribution');
%         title(str);
%         display(t);
%     end
%     if (abs(t-t_end/10) < del_t || abs(t-t_end/5) < del_t || abs(t - t_end) < 2*del_t && false)
%         figure
%         hold on
%         contourf(T_old,20);
%         xlabel('Non- dimensional y-distance');
%         ylabel('Non- dimensional x-distance');
%         zlabel('Non- dimensional temperature');
%         str = sprintf('Temperature profile at time = %f',t);
%         title(str);
%         display(t);
%     end
    %%
    % Calculating the T_max
    %% ADI method
    % In x-direction: [I - k*del_x^2]T` = T^(n)
    % In y-direction: [I - k*del_y^2]T = T`
    % Uses the Thomas algorithm to solve the euler implicit equation as above
    %%
    % Step 1
    % Declaring the matrix A1 and vector b1
    % b is the same as the temperature in vector form
    b1 = mat2vecx(T_old,N_x,N_y);
    A1 = zeros(N_x*N_y,3);
    for n = 1: (N_x *N_y)
        if(mod(n,N_x) == 0)
            % boundary with T=1
            A1(n,2) = 1;
            A1(n,1) = 0;
            A1(n,3) = 0;
            b1(n) = 1;
        elseif(mod(n-1,N_x) == 0 || n <= N_x || n > (N_y-1)*N_x)
            % other boundaries, T=0
            A1(n,2) = 1;
            A1(n,1) = 0;
            A1(n,3) = 0;
            b1(n) = 0;
        else
            A1(n,1) = -K(n)*del_t/del_x^2;
            A1(n,2) = 1+2*K(n)*del_t/del_x^2;
            A1(n,3) = -K(n)*del_t/del_x^2;
        end
    end
    A1(1,1) = 0;
    A1(N_x*N_y, 3) = 0;
    T_star = Thomas_algo(A1,b1,N_x*N_y);
    T_star = vecx2mat(T_star,N_x,N_y);
    
    %%
    % Step 2
    % Declaring the matrix A2 and vector b2
    % b is the same as the intermediate solved temperature in vector form
    b2 = mat2vecy(T_star,N_x,N_y);
    A2 = zeros(N_x*N_y,3);
    for n = 1: (N_x *N_y)
        if(n > (N_x-1) * N_y)
            % boundary with T=1
            A2(n,2) = 1;
            A2(n,1) = 0;
            A2(n,3) = 0;
            b2(n) = 1;
        elseif(mod(n,N_y) == 0 || mod(n-1,N_y) == 0 || n < N_y)
            % boundaries with T=0
            A2(n,2) = 1;
            A2(n,1) = 0;
            A2(n,3) = 0;
            b2(n) = 0;
        else
            % bulk
            A2(n,1) = -K(n)*del_t/del_y^2;
            A2(n,2) = 1+2*K(n)*del_t/del_y^2;
            A2(n,3) = -K(n)*del_t/del_y^2;
        end
    end
    A2(1,1) = 0;
    A2(N_x*N_y, 3) =0;
    T_new = Thomas_algo(A2,b2,N_x*N_y);
    T = vecy2mat(T_new,N_x,N_y);
    
    
    %% Iteration over the same variables
    t = t + del_t;
    iter = iter + 1;
    T_old = T;
    
    %% convergence check
    if(min(min(T_old)) < 0)
        str1 = sprintf('divergence observed at iteration number %d', iter-1);
        display(str1);
        break;
    end
end
display(iter)
T_resize = imresize3(T_solution, [64, 64, 64]);
% figure;
% contourslice(double(T_solution), [], linspace(1, size(T_solution,2), 8), []);
% xlabel('X');
% ylabel('y');
% zlabel('Time');
% view(3)
% axis tight
% title('Contour slices for original solution');
% figure
% contourslice(double(T_resize), [], linspace(1, size(T_resize,2), 8), []);
% xlabel('X');
% ylabel('y');
% zlabel('Time');
% view(3)
% axis tight
% title('Contour slices for interpolated/resized solution')
end
