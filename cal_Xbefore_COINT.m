
function [Xbefore, u_opt, du_opt, u_optstar, motor_cost, energy_cost] = cal_Xbefore_COINT(pr, W, C, xstar, delay, initial_cond_0, init_1, init_2, count, Xbefore, h_pt, pert_t)
% W: matrix, (NN, NN)
% xstar: matrix, (NN, ntrial)
% xstar_k: column vector, (NN, 1)
% Xbefore: cell, (1, ntrial). Xbefore{i} is a matrix, (NN, n_timepoints)

% Calculate parameter P for u(t)
I = eye(size(W));
A = W - I;
Q = lyap(A', C' * C);
B = ones(size(A)) / pr.lambda;
P = are(A, B, Q); % Algebraic Riccati Equation solution

% Calculate Y
Acl = A - P / pr.lambda;
BB = (P * P) / (pr.lambda ^ 2);
Y = lyap(Acl', BB);

% Define time intervals
delta_t = 100; % Target speed = 250 degree/sec
x_inter = pr.tfinal / delta_t;
ahead = pr.t_move / delta_t + 1;

tspan = cell(1, x_inter);
for i = 1:x_inter
    tspan{i} = linspace(delta_t * (i - 1), delta_t * i, delta_t);
end
tspan{x_inter} = linspace(pr.tfinal - delta_t, delay, delay - (pr.tfinal - delta_t));

% Set initial condition
if count == 1
    initial_cond = initial_cond_0;
else
    initial_cond = Xbefore{count - 1}(:, end);
end

% Determine xstar_k and calculate Xbefore
if count <= 2
    xstar_k = xstar(:, init_2);
    Xbefore{count} = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t, 'track');
else
    if init_1 == init_2 % Static
        xstar_k = xstar(:, init_2);
        Xbefore{count} = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t, 'static');
    else % Moving
        init_ahead = mod(init_2 + ahead - 1, 12) + 1;
        xstar_k = xstar(:, init_ahead);
        Xbefore{count} = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t, 'moving');
    end
end

% Transpose Xbefore for correct orientation
Xbefore{count} = Xbefore{count}';

% Calculate motor and energy costs
dx = Xbefore{count} - xstar_k;
motor_cost = zeros(1, size(dx, 2));
energy_cost = zeros(1, size(dx, 2));

for k = 1:size(dx, 2)
    motor_cost(k) = transpose(dx(:, k)) * (P - pr.lambda * Y) * dx(:, k);
    energy_cost(k) = transpose(dx(:, k)) * Y * dx(:, k);
end

% Calculate inputs
h_bar = pr.x_sp - W * max(0, pr.x_sp);
u_opt = xstar_k - W * max(0, xstar_k) - h_bar - (P * (Xbefore{count} - xstar_k)) / pr.lambda;
du_opt = - (P * (Xbefore{count} - xstar_k)) / pr.lambda;
u_optstar = xstar_k - W * max(0, xstar_k);

if sum(h_pt) == 0
    if count > 2 && init_1 == init_2 % Static
        u_naive = xstar_k - W * max(0, xstar_k) - h_bar;
        u_opt = repmat(u_naive, 1, 100);
    end
else
    if count > 2 && count <= pert_t / delta_t + 1 && init_1 == init_2 % Static during perturbation
        u_naive = xstar_k - W * max(0, xstar_k) - h_bar;
        u_opt = repmat(u_naive, 1, 100);
    end
end

% Nested function to integrate dynamics
function X = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan, h_pt, pert_t, mode)
    if strcmp(mode, 'static')
        if sum(h_pt) == 0
            [~, X] = ode45(@(t, X) rate_dynamics_naive(t, X, pr, W, xstar_k), tspan, initial_cond);
        elseif count < pert_t / 100 + 1
            [~, X] = ode45(@(t, X) rate_dynamics_naive(t, X, pr, W, xstar_k), tspan, initial_cond);
        elseif count == pert_t / 100 + 1
            [~, X] = ode45(@(t, X) rate_dynamics_naive_pt(t, X, pr, W, xstar_k, h_pt), tspan, initial_cond);
        else
            [~, X] = ode45(@(t, X) rate_dynamics_lqr(t, X, pr, W, xstar_k, P), tspan, initial_cond);
        end
    else % Moving
        if count == pert_t / 100 + 1
            [~, X] = ode45(@(t, X) rate_dynamics_lqr_pt(t, X, pr, W, xstar_k, P, h_pt), tspan, initial_cond);
        else
            [~, X] = ode45(@(t, X) rate_dynamics_lqr(t, X, pr, W, xstar_k, P), tspan, initial_cond);
        end
    end
end

% ODE functions
function x_dot = rate_dynamics_lqr(t, X, pr, W, xstar_k, P)
    h_bar = pr.x_sp - W * max(0, pr.x_sp);
    x_dot = (1 / pr.tau) * (-X + W * max(0, X) + h_bar + xstar_k - W * max(0, xstar_k) - h_bar - (P * (X - xstar_k)) / pr.lambda);
end

function x_dot = rate_dynamics_naive(t, X, pr, W, xstar_k)
    h_bar = pr.x_sp - W * max(0, pr.x_sp);
    x_dot = (1 / pr.tau) * (-X + W * max(0, X) + h_bar + xstar_k - W * max(0, xstar_k) - h_bar);
end

function x_dot = rate_dynamics_naive_pt(t, X, pr, W, xstar_k, h_pt)
    h_bar = pr.x_sp - W * max(0, pr.x_sp);
    x_dot = (1 / pr.tau) * (-X + W * max(0, X) + h_bar + h_pt + xstar_k - W * max(0, xstar_k) - h_bar);
end

function x_dot = rate_dynamics_lqr_pt(t, X, pr, W, xstar_k, P, h_pt)
    h_bar = pr.x_sp - W * max(0, pr.x_sp);
    x_dot = (1 / pr.tau) * (-X + W * max(0, X) + h_bar + h_pt + xstar_k - W * max(0, xstar_k) - h_bar - (P * (X - xstar_k)) / pr.lambda);
end

end
