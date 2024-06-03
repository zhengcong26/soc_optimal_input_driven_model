
function [c, xstars] = prms_v1(pr, W, ystar, c0, xstars0)
% pr: struct, parameters for the problem
% W: matrix, (NN, NN)
% ystar: cell, target trajectories
% c0: initial condition for c
% xstars0: initial condition for xstars

n_obs = pr.NN;
n_mov = pr.ntrial;
xstars_std = 0.2;
sp = pr.x_sp;

% Initialize xstars0 and c0 if not provided
if nargin < 4
    xstars0 = randn(n_obs, n_mov) * 0.1 / sqrt(n_mov);
    c0 = randn(2, pr.NN) * 0.1 / sqrt(2);
end

xstars0_i = xstars0';
x0 = [c0; xstars0_i];

% Calculate the parameter matrices
I = eye(size(W));
Q = lyap((W - I)', 2 * I);
[eig_vector, ~] = eig(Q); 
a = fliplr(eig_vector(:, end-(n_obs-1):end)); 

% Ocaml parameters
n_muscles = 2;
m1_slice = pr.NN;

lambda_traj = 1 / n_mov;
lambda_reg = 0.001 / (m1_slice * n_muscles);
gamma = I(1:m1_slice, :);

% Define the cost function
func = @(x) calculate_cost(x, pr, W, ystar, a, lambda_reg, lambda_traj, gamma, sp);

% Set optimization options and run optimization
options = optimoptions('fminunc', 'Display', 'iter', 'MaxIter', 500000, 'MaxFunEvals', 500000);
x = fminunc(func, x0, options);

% Extract results
c = x(1:2, :);
xstars = transpose(x(3:end, :));

end

% Calculate the cost function
function cost = calculate_cost(x, pr, W, ystar, a, lambda_reg, lambda_traj, gamma, sp)
    x = unpack(x, a, pr.NN, pr.ntrial, pr.x_sp, xstars_std, gamma);
    
    C = x(1:2,:);
    x_star = transpose(x);
    
    cost_reg = lambda_reg * norm(C, 'fro')^2;
    cost_traj = 0;
    
    for i = 1:length(ystar)
        Xafter = X_after(pr, W, x_star(:, i));        
        cost_traj = cost_traj + lambda_traj * norm(C * max(0, Xafter) - ystar{i}, 'fro')^2;
    end
    
    cost = cost_reg + cost_traj;
end

% Unpack the optimization variables
function x = unpack(x, a, n_obs, n_mov, sp, xstars_std, gamma)
    C = x(1:2, :);
    x_star_raw = transpose(x(3:end, :));
    
    top_obs = a;
    z = n_obs * n_mov * xstars_std^2;
    x_star = top_obs * x_star_raw;
    x_star = sqrt(z / norm(x_star, 'fro')^2) * x_star; 
    
    x_motor = [x_star + sp, sp];
    x_motor = gamma * x_motor;
    h = linsolve(x_motor' * x_motor, x_motor');
    C = C - C * x_motor * h;
    
    x = [C;(x_star+sp)'];
end

% Calculate X_after
function X_aft = X_after(pr, W, xstar)
    tspan = linspace(0, pr.t_move, pr.n_move);
    h = pr.h;
    initial_cond = xstar;

    [~, Xafter] = ode45(@(t, X) rate_dynamics2_ode(t, X, pr, W, h), tspan, initial_cond);
    X_aft = Xafter';
end

% Define the ODE for X_after
function x_dot = rate_dynamics2_ode(t, X, pr, W, h)
    h_bar = pr.x_sp - W * max(0, pr.x_sp);
    x_dot = (1 / pr.tau) * (-X + W * max(0, X) + h_bar + interp1(tspan, h, t));
end
