
function plot_u(u,d,colo,coint)

if coint == 0
    % CO
    ntrial = size(u,2);
    X_plt_pst = 1*ntrial-(ntrial-1) : 1*ntrial;
    
%     figure
    for i = 1%:ntrial
%         subplot(1,ntrial,X_plt_pst(i))
        hold on
        u{i} = [repmat(u{i}(:,1),1,100), u{i}];
        t=linspace(-100,size(u{i},2)-100,size(u{i},2));
        plot(t',u{i}(d,:)','LineWidth',1.5)
%         plot(t',u{i}(1:d,:)','LineWidth',1)
        %         plot(repmat(t',1,3),u_opt{i}(1:3,:)','LineWidth',1)
        
        box off
%         ylim([-300 600]);
        %         xlim([-50 1100]);
        if i == 1
            ylabel('u activity');
        else
            %         xticks([]);
            %         yticks([]);
        end
    end
    
else
    % INT
%     figure
    hold on
    u = [repmat(u(:,1),1,100), u];
    
    t=linspace(-100,size(u,2)-100,size(u,2));
    plot(t',u(d,:)','color',colo,'LineWidth',1)
%     plot(t',u(1:d,:)','LineWidth',1)
    
    box off
%     ylim([-300 600]);
    ylabel('u activity');
    
end

end