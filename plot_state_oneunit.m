

function plot_state_oneunit(pr, X_check,ntrial,num,colo)

hold on

X_check{ntrial} = [repmat(X_check{ntrial}(:,1),1,100), X_check{ntrial}];

t=linspace(-100,size(X_check{ntrial},2)-100,size(X_check{ntrial},2));

plot(t,X_check{ntrial}(num,:)','color',colo,'LineWidth',2) 

box off
ylabel('firing rates (Hz)');
xlim([-150 inf]);
ylim([17 25]);
y1 =  ylim;
xticks([0 100 200 300 400 500 600 700 800 900 1000])
% yticks([15 20 25])
set(gca,'TickDir','out');

h = gca;
% h.XAxis.Visible = 'off';
h.FontSize = 15;
h.LineWidth = 1;
h.TickLength=[0.02,0.025];
% set(gca,'xticklabel',[])
% xticks([]);
% axis off

end