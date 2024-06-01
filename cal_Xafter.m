function Xafter = cal_Xafter(pr, W, xstar)
% W£ºmatrix£¬(NN,NN)
% xstar£ºmatrix£¬(NN,ntrial)
% Xafter£ºcell£¬(1,ntrial)¡£Xafter{i} is a matrix£¬(NN,n_timepoints)

h = pr.h;
% tspan = linspace(0, pr.tfinal, pr.n_timepoints);
tspan = linspace(0, pr.t_move,pr.n_move); % 20230526 constant MT

for i = 1:pr.ntrial
    initial_cond = xstar(:,i);
    Xafter{i} = integrate_dynamics2(pr, W, h, initial_cond);
    Xafter{i} = Xafter{i}'; % transpose, so that the time is in the direction of the ROW vector
end

    % slove ODE 2
    function X = integrate_dynamics2(pr, W, h, initial_cond)
        [~,X] = ode45(@(t,X) rate_dynamics2_ode(t, X, pr, W, h), ...
            tspan, initial_cond);
    end
    
    % define ODE 2£¨after MO£©
    function x_dot = rate_dynamics2_ode(t, X, pr, W, h)
        h_bar = pr.x_sp - W*max(0,pr.x_sp); % phi(x)=max(x,0)
        x_dot = 1/pr.tau* ...
            ( ...
            -X ...
            +W*max(0,X) ...
            +h_bar ...
            +interp1(tspan,h,t) ... % no u(t) after MO
            );
    end

end