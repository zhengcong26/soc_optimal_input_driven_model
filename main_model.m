
%% This code reproduces the work of (Ta-Chu Kao et al., 2021)
%% Written by Cong Zheng

clc
clear

%% step1 optimization of xstar and C, which requires relatively long time
% pr = def_params(); % define parameters
% pr.ntrial = 12; 
% v0 = cal_v0(pr); % velocity parameter of the hand
% hand = cal_hand(pr, v0); % hand trajectory (ground truth, cong)
% pr.lambda = 1;
% load Wsoc21
% %% xstar and C
% % [c,xstars] = prms_v1(pr, Wsoc, hand); % step 1
% [c,xstars] = prms_v1(pr, Wsoc, hand, c, xstars); % step 2-n
% [C,xstar] = Unpack(c,xstars,Wsoc,pr);
% Xafter = cal_Xafter(pr, Wsoc, xstar);
% %% optimize C further
% CC = cal_C(pr, Xafter, xstar, C, hand);
% C = CC;
% %% modify the quadrants
% order = [8,7,6,5,4,3,2,1,12,11,10,9];
% for i = 1:length(order)
%     Hand{i} = hand{order(i)};
%     Xstar(:,i) = xstar(:,order(i));
% end
%% skip step1 and load optimized parameters
load('optimized_xstar_C.mat')

%% step2 perturb 

pr.tfinal = 700;
pr.t_move = 300;
delta_t = 100;

X_check_CO=[];
Y_check_CO=[];
u_CO=[];
motor_cost_CO=[];
energy_cost_CO=[];

X_check_INT=[];
Y_check_INT=[];
u_INT=[];
motor_cost_INT=[];
energy_cost_INT=[];

r2_co_st_rt = [];
rmse_co_st_rt = [];
r2_int_st_rt = [];
rmse_int_st_rt = [];

delay = 700:20:1000;    
pert_amp = [0 20 30 40];

for n = 1:length(pert_amp)
    for init = 1:12
        for rep = 1:10
            
            initial_cond_0 = normrnd(20,1,pr.NN,1);
            
            pert_t = 500;
            pert_index = randsample(200,100);
            h_pt = ones(200,1)*pert_amp(n);
            h_pt(pert_index)=0;
            
            
            for d = 1:length(delay)
                
                Xbefore_opt_CO = cell(1,7);
                Xbefore_opt_INT = cell(1,7);
                motor_cost_opt_CO = [];
                energy_cost_opt_CO = [];
                u_opt_CO = [];
                du_opt_CO = [];
                u_optstar_CO = [];
                motor_cost_opt_INT = [];
                energy_cost_opt_INT = [];
                u_opt_INT = [];
                du_opt_INT = [];
                u_optstar_INT = [];
                X_check_opt_CO = [];
                X_check_opt_INT = [];
                Y_check_opt_CO = [];
                Y_check_opt_INT = [];
                
                tic
                %% static
                for i = 1:pr.tfinal/100
                    init_1 = init;
                    init_2 = init;
                    
                    [Xbefore_opt_CO, u_opt_CO{i}, du_opt_CO{i}, u_optstar_CO{i}, motor_cost_opt_CO{i}, energy_cost_opt_CO{i}] = ...
                        cal_Xbefore_COINT(pr, Wsoc, CC1, Xstar, delay(d), initial_cond_0, init_1, init_2, i, Xbefore_opt_CO, h_pt, pert_t);
                    
                end
                
                init_final = init;
                
                [X_check_opt_CO,Y_check_opt_CO,r2_co_st_rt(n,init,rep,d),rmse_co_st_rt(n,init,rep,d)] = simu_X_INT_plot(pr, Wsoc, CC1, Hand, Xbefore_opt_CO, Xafter, init, init_final);
                X_check_CO(n,init,rep,d,:,:)=X_check_opt_CO;
                Y_check_CO(n,init,rep,d,:,:)=Y_check_opt_CO;
                u_CO(n,init,rep,d,:,:)=cell2mat(u_opt_CO);
                motor_cost_CO(n,init,rep,d,:,:)=cell2mat(motor_cost_opt_CO);
                energy_cost_CO(n,init,rep,d,:,:)=cell2mat(energy_cost_opt_CO);
%             

            close all
            
            %% moving
            for i = 1:pr.tfinal/100
                if i == 1
                    init_1 = init+i-1;
                    if init_1 > 12
                        init_1 = init_1-12;
                    end
                    init_2 = init_1;
                else
                    init_2 = init+i-1;
                    if init_2 > 12
                        init_2 = init_2-12;
                    end
                    init_1 = init_2-1;
                    if init_1 <=0
                        init_1 = 12;
                    end
                end
                
                [Xbefore_opt_INT, u_opt_INT{i}, du_opt_INT{i}, u_optstar_INT{i}, motor_cost_opt_INT{i}, energy_cost_opt_INT{i}] = ...
                    cal_Xbefore_COINT(pr, Wsoc, CC1, Xstar, delay(d), initial_cond_0, init_1, init_2, i, Xbefore_opt_INT, h_pt, pert_t);
                
            end
            
            init_final = init_2+4;
            if init_final > 12
                init_final = init_final-12;
            end
            
            [X_check_opt_INT,Y_check_opt_INT,r2_int_st_rt(n,init,rep,d),rmse_int_st_rt(n,init,rep,d)] = simu_X_INT_plot(pr, Wsoc, CC1, Hand, Xbefore_opt_INT, Xafter, init, init_final);
            X_check_INT(n,init,rep,d,:,:)=X_check_opt_INT;
            Y_check_INT(n,init,rep,d,:,:)=Y_check_opt_INT;
            u_INT(n,init,rep,d,:,:)=cell2mat(u_opt_INT);
            motor_cost_INT(n,init,rep,d,:,:)=cell2mat(motor_cost_opt_INT);
            energy_cost_INT(n,init,rep,d,:,:)=cell2mat(energy_cost_opt_INT);
            
            close all
            toc
            end
        end
    end
end



%% step3 calculate and plot
%% calculate RT

for k = 1:size(r2_co_st_rt,1)
    for i = 1:size(r2_co_st_rt,2)
        for j = 1:size(r2_co_st_rt,3)
            try
                rt_1(k,i,j)=find(squeeze(r2_co_st_rt(k,i,j,:))>0.9,1);
            catch
                rt_1(k,i,j)=size(r2_co_st_rt,4); %
            end
        end
    end
end

for k = 1:size(r2_int_st_rt,1)
    for i = 1:size(r2_int_st_rt,2)
        for j = 1:size(r2_int_st_rt,3)
            try
                rt_2(k,i,j)=find(squeeze(r2_int_st_rt(k,i,j,:))>0.9,1);
            catch
                rt_2(k,i,j)=size(r2_int_st_rt,4); %
            end
        end
    end
end

rt_co=(rt_1-1)*20;
rt_int=(rt_2-1)*20;

% rt_co_1 = (rt_co)'
% rt_int_1 = (rt_int)'

rt_co_1 = (reshape(rt_co, 4, 120))';
rt_int_1 = (reshape(rt_int, 4, 120))';

% s

%% plot hand

target = Y_check_CO;

% hand in one fig
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

figure
hold on
for i = 1:size(target,1)
    y = target(i,:); 
    plot_hand_y(Hand, y, colo(1,:), 1)
end
set(gcf,'position',[30,300,300,300]);

%% plot perturb hand
target = Y_check_INT{2,1};
% hand in one fig
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

figure
plot_hand_y(Hand, target, colo(2,:), 3)
set(gcf,'position',[30,300,300,300]);

%% plot PSTH

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
endColor = colo(2,:);
c1 = gradient_ramp(endColor,pr.ntrial,3);

X_check = X_check_INT(1,:);

figure
for k = 181:200
    subplot(5,4,k-180)
    for i = 1:12
        plot_state_oneunit(pr,X_check,i,k,c1(i,:))
        title(num2str(k))
    end
    %     set(gcf,'position',[30,300,400,250]);
end

% select one
figure
set(gcf,'position',[30,100,350,200]);
k = 197;
for i = 1:12
    plot_state_oneunit(pr,X_check,i,k,c1(i,:))
    %     title(num2str(k))
end
% set(gca,'XColor', 'none')

%% plot u_opt

u = u_INT{1,1};
% u = u_CO{1,1};

it = [102 120 167]; 

c1 = copper(4);

figure
hold on
for i = 1:length(it)
    plot_u(u,it(i),c1(i,:),1)
end
set(gcf,'position',[30,100,230,80]);
set(gca,'yscale','log')
ylim([-100 6000]);
box off
set(gca,'TickDir','out');
h = gca;
h.FontSize = 8.5;
h.TickLength = [0.025,0.03];
xticks([0 100 200 300 400 500 600 700])

%% plot motor_cost
motor_cost_CO_n = vertcat(motor_cost_CO{:});
motor_cost_INT_n = vertcat(motor_cost_INT{:});

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

figure
hold on

% moving
y = (motor_cost_INT_n)';

t = 1:1:size(y,1);
y0 = mean(y,2,'omitnan');
y1 = y0';

for i = 1:20%size(y,2)
    plot(t,y(:,i),'color',[0.80392	0.33333	0.33333 0.2],'LineStyle',linestyle{1},'LineWidth',1); 
end

plot(t,y1,'color',colo(2,:),'LineStyle',linestyle{1},'LineWidth',2); % black

std2 = [];
std0 = std(y,0,2,'omitnan')/sqrt(length(y(1,:))-1); % std
std1 = [y0-std0, y0(end:-1:1)+std0(end:-1:1)];
std2(1,1:length(y1)) = std1(1:length(y1),1)';
std2(1,length(y1)+1:2*length(y1)) = std1(1:length(y1),2)';

if std0 ~= 0
    h=fill([t, t(end:-1:1)],std2,'b');
    set(h,'FaceColor',colo(2,:),'FaceAlpha',0.2,'EdgeColor','none');
end

% static
y = motor_cost_CO_n';

t = 1:1:size(y,1);
y0 = mean(y,2,'omitnan');
y1 = y0';

for i = 1:10%size(y,2)
    plot(t,y(:,i),'color',[0.06275	0.30588	0.5451 0.2],'LineStyle',linestyle{1},'LineWidth',1); 
end

plot(t,y1,'color',colo(1,:),'LineStyle',linestyle{1},'LineWidth',2); 

std2 = [];
std0 = std(y,0,2,'omitnan')/sqrt(length(y(1,:))-1); % std
std1 = [y0-std0, y0(end:-1:1)+std0(end:-1:1)];
std2(1,1:length(y1)) = std1(1:length(y1),1)';
std2(1,length(y1)+1:2*length(y1)) = std1(1:length(y1),2)';

if std0 ~= 0
    h=fill([t, t(end:-1:1)],std2,'b');
    set(h,'FaceColor',colo(1,:),'FaceAlpha',0.2,'EdgeColor','none');
end

% figure setup
set(gca,'yscale','log')
% axis([0 700 -10 100000])
xlabel('Time (ms)')
ylabel('Prospective motor error')
box off
set(gca,'TickDir','out');
h = gca;
h.FontSize = 11;
set(gca,'color','none');
% set(gcf,'position',[30,300,350,210]);
set(gcf,'position',[30,300,1000,500]);
h.TickLength=[0.015,0.025];
h.LineWidth=1;
xlim([-50 700])
% xlim([-50 1000])
ylim([-10 1000000])
% ylim([-10 100000])
xticks([0 100 200 300 400 500 600 700 800 900 1000])
% xticks([0 200 500 1000])
yticks([1 10E4])

%% plot naive and LQR input motor cost
figure
set(gcf,'position',[100,500,200,100]);
plot(motor_cost_CO{1, 1},'color',[0.47843	0.21569	0.5451],'linewidth',1.2)
hold on
plot(motor_cost_CO_naive{1, 1},'color',[0.80392	0.58824	0.80392],'linewidth',1.2)
box off
set(gca,'TickDir','out');
% xlabel('time (ms)')
% ylabel('prospective motor error')
set(gca,'color','none');
set(gca,'yscale','log');
h = gca;
xlim([-50 700])
xticks([0 100 200 300 400 500 600 700])
h.FontSize = 11;
set(gca,'color','none');
h.TickLength=[0.02,0.025];


