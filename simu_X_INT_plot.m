
function [X_check_INT,y_check3_INT,r2,rmse] = simu_X_INT_plot(pr, Wsoc, CC1, hand, Xbefore_INT, Xafter, init, init_ahead)

X_INT = [cell2mat(Xbefore_INT), Xafter{init_ahead}];
%% Xafter_simu
Xbefore_end = Xbefore_INT{end}(:,end); % select last bin X to seed move
Xafter_simu_INT = cal_Xafter_INT(pr, Wsoc, Xbefore_end); % cal x after MO
%% X_check
figure
set(gcf,'position',[100,150,200,650]);
X_check_INT = [cell2mat(Xbefore_INT), Xafter_simu_INT];

subplot(3,1,1)
plot_X_INT(X_INT, X_check_INT) % X differs from X_check in after only, compare move

%% y_check
y_check3_INT = CC1 * max(0,Xafter_simu_INT);

subplot(3,1,2)
plot_y_INT(hand{init_ahead}, y_check3_INT) %

[r2,rmse] = cal_r(hand{init_ahead},y_check3_INT);

%% new hand check3 y-hand
subplot(3,1,3)
colo = [0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
plot_hand_y_INT(hand, y_check3_INT, init, init_ahead, colo(2,:), 0)
