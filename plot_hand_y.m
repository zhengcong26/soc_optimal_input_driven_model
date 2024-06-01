
function plot_hand_y(hand, y_check, colo, onefig)
ntrial = size(hand,2);
if onefig == 0
    hand_plt_pst_1 = 3*ntrial-(ntrial-1) : 3*ntrial;
    hand_plt_pst_2 = 4*ntrial-(ntrial-1) : 4*ntrial;
    
    for trial = 1:ntrial
        subplot(4,ntrial,[hand_plt_pst_1(trial), hand_plt_pst_2(trial)])
        for i = 1:ntrial
            plot(hand{i}(1,:),hand{i}(2,:),'k')
            hold on
        end
        plot(y_check{trial}(1,:), y_check{trial}(2,:),'color','r','linewidth',2)
        %         plot(hand_check{trial}(1,:), hand_check{trial}(2,:),'o','color','r')
        
        box off
        axis equal
        axis off
        xlim([-20 20]);
        ylim([-20 20]);
    end
elseif onefig == 1
    endColor = colo;
    c1 = gradient_ramp(endColor,ntrial,5);
%     c1 = copper(ntrial);
%     figure
    hold on
    for trial = 1:ntrial
        plot(hand{trial}(1,:),hand{trial}(2,:),'k:','linewidth',1)
        plot(y_check{trial}(1,:), y_check{trial}(2,:),'color',c1(trial,:),'linewidth',2)
        
        box off
        axis equal
        axis off
        xlim([-20 20]);
        ylim([-20 20]);
    end
    
elseif onefig == 2
    endColor = colo;
    c1 = gradient_ramp(endColor,ntrial,5);
%     c1 = copper(ntrial);
%     figure
    hold on
    for trial = 1:ntrial
        plot(hand{trial}(1,:),hand{trial}(2,:),'k:','linewidth',1)
        plot(squeeze(y_check(trial,1,:)), squeeze(y_check(trial,2,:)),'color',c1(trial,:),'linewidth',2)
        
        box off
        axis equal
        axis off
        xlim([-20 20]);
        ylim([-20 20]);
    end
    
else
    endColor = colo;
    c1 = gradient_ramp(endColor,ntrial,5);
%     c1 = copper(ntrial);
%     figure
    hold on
    for trial = 1:ntrial
        plot(hand{trial}(1,:),hand{trial}(2,:),'k:','linewidth',1.2)
        
        box off
        axis equal
        axis off
        xlim([-20 20]);
        ylim([-20 20]);
    end 
    s = scatter(hand{1}(1,end),hand{1}(2,end), 200,'go','filled');
    s.MarkerFaceColor = [0.13333	0.5451	0.13333];
    s.MarkerEdgeColor = [0.13333	0.5451	0.13333];
    scatter(hand{3}(1,end),hand{3}(2,end), 200,'go','MarkerEdgeColor',[0.13333	0.5451	0.13333])
    
    plot(y_check(1,:), y_check(2,:),'color',c1(trial,:),'linewidth',3)
    
end

end