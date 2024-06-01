function v0 = cal_v0(pr)

tau_reach = pr.tau_reach;
% tfinal = pr.tfinal;
tfinal = pr.t_move;
d = pr.d_reach;

v0 = 1;
fun_v = @(t) v0*(t/tau_reach).^2.*exp(-0.5*(t/tau_reach).^2);
reach_length = integral(fun_v,0,tfinal);
v0 = d/reach_length;

end