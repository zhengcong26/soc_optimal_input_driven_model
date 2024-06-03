%% This code reproduces the work of (Ta-Chu Kao et al., 2021)
%% Written by Cong Zheng

clc;
clear;

%% Step 1: Optimization of xstar and C

% % Step 1-1 Define parameters and calculate initial values
% pr = def_params(); % Define parameters
% 
% % Calculate initial velocity parameter and hand trajectory
% v0 = cal_v0(pr); % Velocity parameter of the hand
% hand = cal_hand(pr, v0); % Hand trajectory (ground truth)
% 
% % Load predefined weight matrix Wsoc21
% load Wsoc21;
% 
% % Step 1-2: Calculate xstar and C, or load optimized parameters
% % If running for the first time, use prms_v1 to calculate initial c and xstars
% % [c, xstars] = prms_v1(pr, Wsoc, hand); % Step 1
% 
% % If initial c and xstars are available, pass them in to continue optimization
% [c, xstars] = prms_v1(pr, Wsoc, hand, c, xstars); % Step 2-n
% 
% % Unpack c and xstars to obtain C and xstar
% [C, xstar] = Unpack(c, xstars, Wsoc, pr);
% 
% % Calculate Xafter
% Xafter = cal_Xafter(pr, Wsoc, xstar);
% 
% % Step 1-3: Further optimize C
% CC = cal_C(pr, Xafter, xstar, C, hand);
% C = CC;
% 
% % Step 4: Adjust the order of quadrants
% order = [8, 7, 6, 5, 4, 3, 2, 1, 12, 11, 10, 9];
% Hand = cell(1, length(order)); % Pre-allocate Hand
% Xstar = zeros(size(xstar, 1), length(order)); % Pre-allocate Xstar
% 
% for i = 1:length(order)
%     Hand{i} = hand{order(i)};
%     Xstar(:, i) = xstar(:, order(i));
% end

%% Step 1: Load optimized parameters directly to skip lengthy optimization process
load('optimized_xstar_C.mat')

%% Step 2: Perturbation simulation setup

% Simulation parameters
pr.tfinal = 700;
pr.t_move = 300;
delta_t = 100;

% Initialize variables for storing results
r2_static = [];
rmse_static = [];
X_check_static = [];
Y_check_static = [];
u_static = [];
motor_cost_static = [];
energy_cost_static = [];
r2_moving = [];
rmse_moving = [];
X_check_moving = [];
Y_check_moving = [];
u_moving = [];
motor_cost_moving = [];
energy_cost_moving = [];

% Perturbation parameters
delay = 700:20:1000;    
pert_amp = [0, 20, 30, 40];
pert_t = 500;

% Main simulation loop
for n = 1:length(pert_amp)
    for init = 1:12
        for rep = 1:10
            
            initial_cond_0 = normrnd(20,1,pr.NN,1);
            pert_index = randsample(200,100);
            h_pt = ones(200,1) * pert_amp(n);
            h_pt(pert_index) = 0;
            
            for d = 1:length(delay)
                %% Static Condition
                [X_check_opt, Y_check_opt, u_opt, motor_cost_opt, energy_cost_opt, ...
                    r2, rmse] = run_simulation(pr, Wsoc, CC1, Xstar, delay(d), ...
                    initial_cond_0, init, pert_t, h_pt, 'static');
                
                r2_static(n, init, rep, d) = r2;
                rmse_static(n, init, rep, d) = rmse;
                X_check_static(n, init, rep, d, :, :) = X_check_opt;
                Y_check_static(n, init, rep, d, :, :) = Y_check_opt;
                u_static(n, init, rep, d, :, :) = cell2mat(u_opt);
                motor_cost_static(n, init, rep, d, :, :) = cell2mat(motor_cost_opt);
                energy_cost_static(n, init, rep, d, :, :) = cell2mat(energy_cost_opt);
                
                %% Moving Condition
                [X_check_opt, Y_check_opt, u_opt, motor_cost_opt, energy_cost_opt, ...
                    r2, rmse] = run_simulation(pr, Wsoc, CC1, Xstar, delay(d), ...
                    initial_cond_0, init, pert_t, h_pt, 'moving');
                
                r2_moving(n, init, rep, d) = r2;
                rmse_moving(n, init, rep, d) = rmse;
                X_check_moving(n, init, rep, d, :, :) = X_check_opt;
                Y_check_moving(n, init, rep, d, :, :) = Y_check_opt;
                u_moving(n, init, rep, d, :, :) = cell2mat(u_opt);
                motor_cost_moving(n, init, rep, d, :, :) = cell2mat(motor_cost_opt);
                energy_cost_moving(n, init, rep, d, :, :) = cell2mat(energy_cost_opt);
            end
        end
    end
end

%% Step 3: Calculate and plot results

% Calculate reaction time (RT)
rt_static = calculate_reaction_time(r2_static);
rt_moving = calculate_reaction_time(r2_moving);

% rt_static_1 = reshape(rt_static, 4, 120)';
% rt_static_1 = reshape(rt_moving, 4, 120)';

% Plot hand trajectories
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
plot_hand_trajectories(Hand, Y_check_static, colo);

%% Function Definitions

function [X_check_opt, Y_check_opt, u_opt, motor_cost_opt, energy_cost_opt, r2, rmse] = run_simulation(pr, Wsoc, CC1, Xstar, delay, ...
    initial_cond_0, init, pert_t, h_pt, mode)
    % Function to run the simulation for both static and moving conditions
    
    Xbefore_opt = cell(1, 7);
    
    for i = 1:pr.tfinal/100
        if strcmp(mode, 'static')
            init_1 = init;
            init_2 = init;
        else % mode is 'moving'
            if i == 1
                init_1 = mod(init + i - 1 - 1, 12) + 1;
                init_2 = init_1;
            else
                init_2 = mod(init + i - 1 - 1, 12) + 1;
                init_1 = mod(init_2 - 2, 12) + 1;
            end
        end
        
        [Xbefore_opt, u_opt{i}, du_opt{i}, u_optstar{i}, motor_cost_opt{i}, energy_cost_opt{i}] = ...
            cal_Xbefore_COINT(pr, Wsoc, CC1, Xstar, delay, initial_cond_0, init_1, init_2, i, Xbefore_opt, h_pt, pert_t);
    end
    
    if strcmp(mode, 'static')
        init_final = init;
    else % mode is 'moving'
        init_final = mod(init_2 + 3, 12) + 1;
    end

    [X_check_opt, Y_check_opt, r2, rmse] = ...
        simu_X_INT_plot(pr, Wsoc, CC1, Hand, Xbefore_opt, Xafter, init, init_final);
end

function rt = calculate_reaction_time(r2)
    % Function to calculate reaction time (RT)
    for k = 1:size(r2, 1)
        for i = 1:size(r2, 2)
            for j = 1:size(r2, 3)
                try
                    rt(k, i, j) = find(squeeze(r2(k, i, j, :)) > 0.9, 1);
                catch
                    rt(k, i, j) = size(r2, 4); %
                end
            end
        end
    end
    rt = (rt - 1) * 20;
end

function plot_hand_trajectories(Hand, target, colo)
    % Function to plot hand trajectories

    figure
    hold on;
    for i = 1:size(target, 1)
        y = target(i, :); 
        plot_hand_y(Hand, y, colo(1, :), 1);
    end
    set(gcf, 'position', [30, 300, 300, 300]);
end
