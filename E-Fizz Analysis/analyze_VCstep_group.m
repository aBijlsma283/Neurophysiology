function analyze_VCstep_group(filename,group1,group2,mode)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_VCstep_group is a function used for a voltage clamp sawtooth protocol
% analysis; (magneto version)

% This function retrives all the files pretagged with "analyzed_" and
% performs group analysis on the variables calculated by the analyze_ST functions. 

close all;

%% import the data from the current folder 
if nargin == 1
    filesgroup1 = dir ('analyzed_VCstep*.mat');
elseif nargin == 4
    filesgroup1 = dir (['analyzed_VCstep*' '_' (group1) '*.mat']);
    filesgroup2 = dir (['analyzed_VCstep*' '_' (group2) '*.mat']);
end

%% Gathering of data
if nargin == 1
    for ctr = 1:size(filesgroup1,1)
        load (filesgroup1(ctr).name, 'output')
        gdata.preprom(:,1) = output.data.plots.pks(:,1);
        gdata.preprom(gdata.preprom(:,1) == 0) = NaN;
        gdata.prom(ctr,1) = nanmean(gdata.preprom(:,1));
        
        gdata.IVgroup1(ctr,:) = gdata.preprom(:,1);
        gdata.stimgroup1(ctr,:) = round(output.Vmon.raw(:,3000));
        
        output.data.plots.pksindex(isnan(output.data.plots.pksindex(:,1)),1) = 1001;
        gdata.pretiming(:,1) = output.data.time(output.data.plots.pksindex(:,1)-1000);
        gdata.pretiming(gdata.pretiming(:,1) == 0) = NaN;
        gdata.timing(ctr,1) = nanmean(gdata.pretiming(:,1));
        
        gdata.events(ctr,1) = mean(output.data.plots.events(:,1));
        
        gdata.preperdiffgroup1(:,1) = output.data.plots.pksindex(:,2)-output.data.plots.pksindex(:,1);
        gdata.preperdiffgroup1(gdata.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.preperdiffgroup1(:,1)))) == 0
            gdata.perdiffgroup1(ctr,1) = output.data.time(round(nanmean(gdata.preperdiffgroup1(:,1))));
        else
            gdata.perdiffgroup1(ctr,1) = NaN;
        end
    end
elseif nargin == 4
    for ctr = 1:size(filesgroup1,1)
        load (filesgroup1(ctr).name, 'output')
        gdata.prom(ctr,1) = nanmean(output.data.plots.pks(:,1));
        gdata.IVgroup1(ctr,:) = output.data.plots.pks(:,1);
        gdata.stimgroup1(ctr,:) = round(output.Vmon.raw(:,output.data.plots.protocol.stimbeginloc + 500));
        output.data.plots.pksindex(isnan(output.data.plots.pksindex(:,1)),1) = output.data.plots.protocol.stimbeginloc + 1;
        gdata.pretiming(:,1) = output.data.time(output.data.plots.pksindex(:,1)-output.data.plots.protocol.stimbeginloc);
        gdata.pretiming(gdata.pretiming(:,1) == 0) = NaN;
        gdata.timing(ctr,1) = nanmean(gdata.pretiming(:,1));
        gdata.latency1(ctr,:) = gdata.pretiming(:,1);
        gdata.events(ctr,1) = mean(output.data.plots.events(:,1));
        gdata.preperdiffgroup1(:,1) = output.data.plots.pksindex(:,2)-output.data.plots.pksindex(:,1);
        gdata.preperdiffgroup1(gdata.preperdiffgroup1(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.preperdiffgroup1(:,1)))) == 0
            gdata.perdiff(ctr,1) = output.data.time(round(nanmean(gdata.preperdiffgroup1(:,1))));
        else
            gdata.perdiff(ctr,1) = NaN;
        end
    end
    
    for ctr = 1:size(filesgroup2,1)
        load (filesgroup2(ctr).name, 'output')
        gdata.prom(ctr,2) = nanmean(output.data.plots.pks(:,1));
        gdata.IVgroup2(ctr,:) = output.data.plots.pks(:,1);
        gdata.stimgroup2(ctr,:) = round(output.Vmon.raw(:,output.data.plots.protocol.stimbeginloc + 500));
        output.data.plots.pksindex(isnan(output.data.plots.pksindex(:,1)),1) = 1001;
        gdata.pretiming(:,1) = output.data.time(output.data.plots.pksindex(:,1)-1000);
        gdata.pretiming(gdata.pretiming(:,1) == 0) = NaN;
        gdata.timing(ctr,2) = nanmean(gdata.pretiming(:,1));
        gdata.latency2(ctr,:) = gdata.pretiming(:,1);
        gdata.events(ctr,2) = mean(output.data.plots.events(:,1));
        gdata.preperdiffgroup2(:,1) = output.data.plots.pksindex(:,2)-output.data.plots.pksindex(:,1);
        gdata.preperdiffgroup2(gdata.preperdiffgroup2(:,1) == 1) = NaN;
        if isnan(round(nanmean(gdata.preperdiffgroup2(:,1)))) == 0
            gdata.perdiff(ctr,2) = output.data.time(round(nanmean(gdata.preperdiffgroup2(:,1))));
        else
            gdata.perdiff(ctr,2) = NaN;
        end
    end
end

gdata.prom(gdata.prom == 0) = NaN;
gdata.timing(gdata.timing == 0) = NaN;
gdata.events(gdata.events == 0) = NaN;
gdata.perdiff(gdata.perdiff == 0) = NaN;

%% boxplots (figure 1)
if nargin == 1
    [gdata.hypo(1),gdata.p(1)] = ttest(gdata.prom(:,1));
    [gdata.hypo(2),gdata.p(2)] = ttest(gdata.timing(:,1));
    [gdata.hypo(3),gdata.p(3)] = ttest(gdata.events(:,1));
    [gdata.hypo(4),gdata.p(4)] = ttest(gdata.perdiffgroup1(:,1));
    
    line_width = 1.5;
    marker_size = 10;
    condition = {filename};
    
    figure(1);
    set(gcf,'position', get(0,'screensize'));
    
    for ctr = 1 : size(filesgroup1,1)
        h1 = subplot(2,2,1); hold on;
        plot(gdata.prom(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Amplitude. p = ' num2str(gdata.p(1))],'FontSize',12); xlabel('Condition'); ylabel('Amplitude [pA]');
        
        h2 = subplot(2,2,2); hold on;
        plot(gdata.timing(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Latency. p = ' num2str(gdata.p(2))],'FontSize',12); xlabel('Condition'); ylabel('Time [s]');
        A = gca;
        A.YAxis.Exponent = 0;
        
        h3 = subplot(2,2,3); hold on;
        plot(gdata.events(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['Number of Events. p = ' num2str(gdata.p(3))],'FontSize',12); xlabel('Condition'); ylabel('Events');
        
        h4 = subplot(2,2,4); hold on;
        plot(gdata.perdiffgroup1(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1-2 Peak Latency Difference. p = ' num2str(gdata.p(4))],'FontSize',12); xlabel('Condition'); ylabel('Time [s]');
    end
    
    boxplot_position = [ones(size(gdata.prom,1), 1)*0.5, ones(size(gdata.prom,1), 1)*2.5];
    boxplot_width = [0.3, 0.3];
    
    axes(h1); hold on
    boxplot(gdata.prom(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    xlim([0 1.5]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h2); hold on
    boxplot(gdata.timing(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    xlim([0 1.5]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h3); hold on
    boxplot(gdata.events(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    xlim([0 1.5]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h4); hold on
    boxplot(gdata.perdiffgroup1(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    xlim([0 1.5]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
elseif nargin == 4

    if strcmp(mode, 'Yes') 
        [gdata.hypo(1),gdata.p(1)] = ttest(gdata.prom(:,1),gdata.prom(:,2));
        [gdata.hypo(2),gdata.p(2)] = ttest(gdata.timing(:,1),gdata.timing(:,2));
        [gdata.hypo(3),gdata.p(3)] = ttest(gdata.events(:,1),gdata.events(:,2));
        [gdata.hypo(4),gdata.p(4)] = ttest(gdata.perdiff(:,1),gdata.perdiff(:,2));
        plotstyle = '-ob';
    else
        [gdata.hypo(1),gdata.p(1)] = ttest2(gdata.prom(:,1),gdata.prom(:,2));
        [gdata.hypo(2),gdata.p(2)] = ttest2(gdata.timing(:,1),gdata.timing(:,2));
        [gdata.hypo(3),gdata.p(3)] = ttest2(gdata.events(:,1),gdata.events(:,2));
        [gdata.hypo(4),gdata.p(4)] = ttest2(gdata.perdiff(:,1),gdata.perdiff(:,2));
        plotstyle = 'ob';
    end
    
    line_width = 1.5;
    marker_size = 10;
    condition = {group1, group2};
    
    figure(1);
    set(gcf,'position', get(0,'screensize'));
    
    for ctr = 1 : size(filesgroup1,1)
        h1 = subplot(2,2,1); hold on;
        plot(gdata.prom(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Amplitude. p = ' num2str(gdata.p(1))],'FontSize',12); xlabel('Condition'); ylabel('Amplitude [pA]');
        
        h2 = subplot(2,2,2); hold on;
        plot(gdata.timing(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Latency. p = ' num2str(gdata.p(2))],'FontSize',12); xlabel('Condition'); ylabel('Time [s]');
        A = gca;
        A.YAxis.Exponent = 0;
        
        h3 = subplot(2,2,3); hold on;
        plot(gdata.events(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['Number of Events. p = ' num2str(gdata.p(3))],'FontSize',12); xlabel('Condition'); ylabel('Events');
        
        h4 = subplot(2,2,4); hold on;
        plot(gdata.perdiff(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1-2 Peak Latency Difference. p = ' num2str(gdata.p(4))],'FontSize',12); xlabel('Condition'); ylabel('Time [s]');
    end
    
    boxplot_position = [ones(size(gdata.prom,1), 1)*0.5, ones(size(gdata.prom,1), 1)*2.5];
    boxplot_width = [0.3, 0.3];
    
    axes(h1); hold on
    boxplot(gdata.prom(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.prom(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h2); hold on
    boxplot(gdata.timing(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.timing(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h3); hold on
    boxplot(gdata.events(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.events(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
    
    axes(h4); hold on
    boxplot(gdata.perdiff(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.perdiff(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition)
end 

set(gcf,'Position', get(0, 'Screensize'));
savefig([filename '_Boxplots_figure1']);
saveas(gcf,[filename '_Boxplots_figure1.png']);

%% I/V Curve (fig 2)

if nargin == 4
    gdata.IVgroup1(isnan(gdata.IVgroup1)) = 0;
    gdata.IVgroup2(isnan(gdata.IVgroup2)) = 0;
    
    gdata.Istd1 = std(gdata.IVgroup1);
    gdata.Istd2 = std(gdata.IVgroup2);
    
    gdata.Imean1 = mean(gdata.IVgroup1);
    gdata.Imean2 = mean(gdata.IVgroup2);
    
    gdata.Idifconditions = abs(gdata.Imean1-gdata.Imean2);
    
    figure(2)
    hold on; grid on
    plot(gdata.stimgroup1(1,:),-gdata.Imean1,'r.');
    plot(gdata.stimgroup2(1,:),-gdata.Imean2,'k.');
    plot(gdata.stimgroup1(1,:),gdata.Idifconditions,'b.');
    f = fit(gdata.stimgroup1(1,:)',-gdata.Imean1','smoothingspline','SmoothingParam',0.1);
    r1 = plot(f,'r');
    hold on
    f = fit(gdata.stimgroup2(1,:)',-gdata.Imean2','smoothingspline','SmoothingParam',0.1);
    r2 = plot(f,'k');
    hold on
    f = fit(gdata.stimgroup1(1,:)',gdata.Idifconditions','smoothingspline','SmoothingParam',0.1);
    r3 = plot(f,'b');
    legend off
    xlim([min(gdata.stimgroup1(1,:)) max(gdata.stimgroup1(1,:))])
    
    errorbar(gdata.stimgroup1(1,:),-gdata.Imean1,gdata.Istd1,'r.')
    errorbar(gdata.stimgroup2(1,:),-gdata.Imean2,gdata.Istd2,'k.')
    
    xlabel('Hold potential [mV]','FontSize',16)
    ylabel('Current [nA]','FontSize',16)
    title('I/V curve','FontSize',16)
    legend([r1 r2 r3],{group1,group2,'Residuals'})
    
    set(gcf,'Position', get(0, 'Screensize'));
    savefig([filename '_IVcurve_figure2']);
    saveas(gcf,[filename '_IVcurve_figure2.png']);
    
% L/V Curve (fig 3)

    gdata.latency1(isnan(gdata.latency1)) = 0;
    gdata.latency2(isnan(gdata.latency2)) = 0;
    
    gdata.Lstd1 = std(gdata.latency1);
    gdata.Lstd2 = std(gdata.latency2);
    
    gdata.Lmean1 = mean(gdata.latency1);
    gdata.Lmean2 = mean(gdata.latency2);
    
    gdata.Ldifconditions = abs(gdata.Lmean1-gdata.Lmean2);
    figure(3)
    hold on; grid on
    plot(gdata.stimgroup1(1,:),gdata.Lmean1,'r.');
    plot(gdata.stimgroup2(1,:),gdata.Lmean2,'k.');
    plot(gdata.stimgroup1(1,:),gdata.Ldifconditions,'b.');
    
    f = fit(gdata.stimgroup1(1,:)',gdata.Lmean1','smoothingspline','SmoothingParam',0.1);
    r1 = plot(f,'r');
    hold on
    f = fit(gdata.stimgroup2(1,:)',gdata.Lmean2','smoothingspline','SmoothingParam',0.1);
    r2 = plot(f,'k');
    hold on
    f = fit(gdata.stimgroup1(1,:)',gdata.Ldifconditions','smoothingspline','SmoothingParam',0.1);
    r3 = plot(f,'b');
    legend off
    xlim([min(gdata.stimgroup1(1,:)) max(gdata.stimgroup1(1,:))])
    
    errorbar(gdata.stimgroup1(1,:),gdata.Lmean1,gdata.Lstd1,'r.')
    errorbar(gdata.stimgroup2(1,:),gdata.Lmean2,gdata.Lstd2,'k.')
      
    xlabel('Hold potential [mV]','FontSize',16)
    ylabel('time [s]','FontSize',16)
    A = gca;
    A.YAxis.Exponent = 0;
    title('L/V curve','FontSize',16)
    legend([r1 r2 r3],{group1,group2,'Residuals'})
    
    save(['completegroup_VCstep_' filename], 'gdata', 'output')
    set(gcf,'Position', get(0, 'Screensize'));
    savefig([filename '_LVcurve_figure3']);
    saveas(gcf,[filename '_LVcurve_figure3.png']);
end

close all
% %% Plotting (figure 2)
% 
% figure();
% set(gcf,'position', get(0,'screensize'));
% 
% load (filesgroup1(4).name, 'output')
% StimulusRange = round(output.Vmon.raw(:,1025));
% mygray = analyze_VCstep_gradientmaker([0.8 0.8 0.8],[0 0 0],20);
% for ctr = 1:20
%     subplot(2,3,1:3);
%     color = mygray;
%     r1 = plot(output.data.time,output.Imon.raw(ctr,:),'color',color(ctr,:));
%     hold on; grid on
% end
% 
% load (filesgroup2(4).name, 'output')
% myred = analyze_VCstep_gradientmaker([1 0.8 0.8],[1 0 0],20);
% for ctr = 1:20
%     subplot(2,3,1:3);
%     hold on
%     color = myred;
%     r2 = plot(output.data.time,output.Imon.raw(ctr,:),'color',color(ctr,:));
% end
% 
% ylabel('Current [nA]','FontSize',12)
% xlabel('Time [s]','FontSize',12)
% title('Representation VC-Step Protocol','FontSize',22)
% legend([r1 r2],{'Magnet OFF','Magnet ON'})
% 
% boxplot_position = [ones(size(gdata.prom,1), 1)*0.5, ones(size(gdata.prom,1), 1)*2.5];
% boxplot_width = [0.3, 0.3];
% 
% axes(h2); hold on
% boxplot(gdata.prom(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
% boxplot(gdata.prom(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
% xlim([0 3]); ylim('auto')
% set(gca,'xtick',[1 2],'xticklabel',condition) 
% 
% axes(h3); hold on
% boxplot(gdata.timing(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
% boxplot(gdata.timing(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
% xlim([0 3]); ylim('auto')
% set(gca,'xtick',[1 2],'xticklabel',condition)  
%  
% axes(h4); hold on
% boxplot(gdata.perdiffWithout(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
% boxplot(gdata.perdiffWith(:,1) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
% xlim([0 3]); ylim('auto')
% set(gca,'xtick',[1 2],'xticklabel',condition) 
% 
% save(['completegroup_VCstep_' filename], 'gdata', 'output')
% set(gcf,'Position', get(0, 'Screensize'));
% savefig([filename '_2']);
