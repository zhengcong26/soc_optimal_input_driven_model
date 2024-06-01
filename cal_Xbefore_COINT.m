
function [Xbefore, u_opt, du_opt, u_optstar, motor_cost, energy_cost] = cal_Xbefore_COINT(pr, W, C, xstar, delay, initial_cond_0, init_1, init_2, count, Xbefore, h_pt, pert_t)
% W£ºmatrix£¬(NN,NN)
% xstar£ºmatrix£¬(NN,ntrial)
% xstar_k£ºvolumn vector£¬(NN,1)
% Xbefore£ºcell£¬(1,ntrial)¡£Xbefore{i} is a matrix£¬(NN,n_timepoints)

% calculate parameter P of u(t)
% A'*P + P*A - lambda-1*P*P+Q = 0
% A'*X + X*A - X*B*X + C = 0
I = eye(size(W));
A = W-I;
Q = lyap(A',C'*C);
B = 1/pr.lambda * ones(size(A)); % ones
P = are(A, B, Q); % Algebraic Riccati Equation solution.

% cal Y
Acl = A-P/pr.lambda;
BB = 1/(pr.lambda)^2*P*P;
Y = lyap(Acl',BB);

%% calculate Xbefore

delta_t = 100; 
x_inter = pr.tfinal/delta_t;
ahead = pr.t_move/delta_t+1;

for i = 1:x_inter
    tspan{i} = linspace(delta_t*(i-1), delta_t*i, delta_t);
end
tspan{x_inter}=linspace(pr.tfinal-delta_t, delay, delay-(pr.tfinal-delta_t));

%% initial condition, xstar
if count == 1
    initial_cond = initial_cond_0;   
else
    initial_cond = Xbefore{count-1}(:,end);
end

if count <= 2
    xstar_k = xstar(:,init_2);
    Xbefore{count} = integrate_dynamics2(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t);
else
    if init_1 == init_2 % static
        xstar_k = xstar(:,init_2);
        
        Xbefore{count} = integrate_dynamics1(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t);
          
    else % moving
        init_ahead = init_2+ahead;
        if init_ahead > 12
            init_ahead = init_ahead-12;
        end
        xstar_k = xstar(:,init_ahead);
        
        Xbefore{count} = integrate_dynamics2(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt, pert_t);

    end
end

%%
Xbefore{count} = Xbefore{count}';

dx = Xbefore{count}-xstar_k;
for k = 1:size(dx,2)
    motor_cost(k) = transpose(dx(:,k))*(P-pr.lambda*Y)*dx(:,k);
    energy_cost(k) = transpose(dx(:,k))*Y*dx(:,k);
end

% cal input
h_bar = pr.x_sp - W*max(0,pr.x_sp); 
u_opt = xstar_k - W*max(0,xstar_k) - h_bar - 1/pr.lambda*P*(Xbefore{count}-xstar_k);  
du_opt = - 1/pr.lambda*P*(Xbefore{count}-xstar_k);
u_optstar = xstar_k - W*max(0,xstar_k);

if sum(h_pt) == 0
    if count > 2 && init_1 == init_2
        u_naive = xstar_k - W*max(0,xstar_k) - h_bar;
        u_opt = u_naive*ones(1,100);
    end
else
    if count > 2 && count <= pert_t/delta_t+1 && init_1 == init_2 
        u_naive = xstar_k - W*max(0,xstar_k) - h_bar;
        u_opt = u_naive*ones(1,100);
    end
end

%% CW
% solve ODE 1
    function X = integrate_dynamics1(pr, W, xstar_k, P, initial_cond, tspan, h_pt, pert_t)
        if sum(h_pt) == 0
            [~,X] = ode45(@(t,X) rate_dynamics_naive_ode(t, X, pr, W, xstar_k, P), ...
                tspan, initial_cond);
        else
            if count < pert_t/100+1
                [~,X] = ode45(@(t,X) rate_dynamics_naive_ode(t, X, pr, W, xstar_k, P), ...
                    tspan, initial_cond);
            elseif count == pert_t/100+1
                [~,X] = ode45(@(t,X) rate_dynamics_naive_pt_ode(t, X, pr, W, xstar_k, P, h_pt), ...
                    tspan, initial_cond);
            else
                [~,X] = ode45(@(t,X) rate_dynamics_lqr_ode(t, X, pr, W, xstar_k, P), ...
                    tspan, initial_cond);
            end
        end
    end

% solve ODE 2
    function X = integrate_dynamics2(pr, W, xstar_k, P, initial_cond, tspan, h_pt, pert_t)
        if count == pert_t/100+1
            [~,X] = ode45(@(t,X) rate_dynamics_lqr_pt_ode(t, X, pr, W, xstar_k, P, h_pt), ...
                tspan, initial_cond);
        else
             [~,X] = ode45(@(t,X) rate_dynamics_lqr_ode(t, X, pr, W, xstar_k, P), ...
                tspan, initial_cond);
        end
    end

% define ODE 2£¨before MO£©
    function x_dot = rate_dynamics_lqr_ode(t, X, pr, W, xstar_k, P)
        h_bar = pr.x_sp - W*max(0,pr.x_sp);
        x_dot = 1/pr.tau* ...
            ( ...
            -X ...
            +W*max(0,X) ...
            +h_bar ...
            +xstar_k - W*max(0,xstar_k) - h_bar ... 
            -1/pr.lambda*P*(X-xstar_k) ...
            );
    end

%     % for naive calculation
    function x_dot = rate_dynamics_naive_ode(t, X, pr, W, xstar_k, P)
        h_bar = pr.x_sp - W*max(0,pr.x_sp);
        x_dot = 1/pr.tau* ...
            ( ...
            -X ...
            +W*max(0,X) ...
            +h_bar ...
            +xstar_k - W*max(0,xstar_k) - h_bar ... 
             );
    end

    function x_dot = rate_dynamics_naive_pt_ode(t, X, pr, W, xstar_k, P, h_pt)
        h_bar = pr.x_sp - W*max(0,pr.x_sp);
        x_dot = 1/pr.tau* ...
            ( ...
            -X ...
            +W*max(0,X) ...
            +h_bar ...
            +h_pt ...              
            +xstar_k - W*max(0,xstar_k) - h_bar ... 
             );
    end

    % ODE pert
    function x_dot = rate_dynamics_lqr_pt_ode(t, X, pr, W, xstar_k, P, h_pt)
        h_bar = pr.x_sp - W*max(0,pr.x_sp);
        x_dot = 1/pr.tau* ...
            ( ...
            -X ...
            +W*max(0,X) ...
            +h_bar ...
            +h_pt ...
            +xstar_k - W*max(0,xstar_k) - h_bar ...
             -1/pr.lambda*P*(X-xstar_k) ...
             );
    end


end




