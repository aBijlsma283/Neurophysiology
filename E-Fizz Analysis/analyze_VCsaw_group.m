function analyze_VCsaw_group(filename,group1,group2,mode)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_VCsaw_group is a function used for a voltage clamp sawtooth protocol
% analysis;

% This function retrives all the files pretagged with "analyzed_" and
% performs group analysis on the variables calculated by the analyze_ST functions.

close all;

%% import the data from the current folder
if nargin == 1
    filesgroup1 = dir ('analyzed_VCsaw*.mat');
elseif nargin == 4
    filesgroup1 = dir (['analyzed_VCsaw*' '_' (group1) '*.mat']);
    filesgroup2 = dir (['analyzed_VCsaw*' '_' (group2) '*.mat']);
end

%% Gathering of data
if nargin == 1
    for ctr = 1:size(filesgroup1,1)
        try
            load (filesgroup1(ctr).name, 'output')
            gdata.prom(ctr,:) = output.data.sweep{1}.findpeaks.sweepprom(:,1);
            gdata.pksinitiationvoltage(ctr,:) = output.Vmon.raw(1,output.data.sweep{1}.plots.pksinitiationindex(1,:));
            gdata.timing(ctr,:) = output.data.sweep{1}.plots.locAfterStim(1,:);
            gdata.perdiff(ctr,:) = output.data.sweep{1}.findpeaks.sweepperdiff(1,:);
            gdata.events(ctr,:) = output.data.sweep{1}.findpeaks.sweepevents(1,:);
            gdata.halfwidth(ctr,:) = output.data.sweep{1}.plots.halfwidth(1,:);
            gdata.fnames1{ctr,1} = filesgroup1(ctr).name;
            gdata.Imon.raw1{ctr} = output.Imon.raw(1,:);
            gdata.time1{ctr} = output.data.sweep{1}.time{1};
        catch
            load (filesgroup1(ctr).name, 'output')
            gdata.prom(ctr,:) = output.data.findpeaks.event1prom(1,:);
            gdata.pksinitiationvoltage(ctr,:) = output.Vmon.raw(1,output.data.plots.pksinitiationindex(1,:));
            gdata.timing(ctr,:) = output.data.plots.locAfterStim(1,:);
            gdata.perdiff(ctr,:) = output.data.findpeaks.perdiff(1,:);
        end
    end
    
elseif nargin == 4
    for ctr = 1:size(filesgroup1,1)
        load (filesgroup1(ctr).name, 'output')
        gdata.prom(ctr,1) = output.data.sweep{1}.findpeaks.sweepprom(1,1);
        gdata.adaptation1(ctr,:) = output.data.sweep{1}.findpeaks.sweepprom(:,1);
        gdata.pksinitiationvoltage(ctr,1) = output.Vmon.raw(1,output.data.sweep{1}.plots.pksinitiationindex(1,1));
        gdata.timing(ctr,1) = output.data.sweep{1}.plots.locAfterStim(1,1);
        gdata.perdiff(ctr,1) = output.data.sweep{1}.findpeaks.sweepperdiff(1,1);
        gdata.events(ctr,1) = output.data.sweep{1}.findpeaks.sweepevents(1,1);
        gdata.halfwidth(ctr,1) = output.data.sweep{1}.plots.halfwidth(1,1);
        gdata.fnames1{ctr,1} = filesgroup1(ctr).name;
        gdata.Imon.raw1{ctr} = output.Imon.raw(1,:);
        gdata.time1{ctr} = output.data.sweep{1}.time{1};
        gdata.integral(ctr,1) = output.data.sweep{1}.plots.integral(1,1);
        gdata.onsettime(ctr,1) = output.data.sweep{1}.plots.onsettime(1,1);
        gdata.offsettime(ctr,1) = output.data.sweep{1}.plots.offsettime(1,1);
        gdata.onsetcurrent(ctr,1) = output.data.sweep{1}.plots.onsetcurrent(1,1);
        gdata.offsetcurrent(ctr,1) = output.data.sweep{1}.plots.offsetcurrent(1,1);
        gdata.onsetslopecurrent(ctr,1) = output.data.sweep{1}.plots.onsetslopecurrent(1,1);
        gdata.offsetslopecurrent(ctr,1) = output.data.sweep{1}.plots.offsetslopecurrent(1,1);
    end
    for ctr = 1:size(filesgroup2,1)
        load (filesgroup2(ctr).name, 'output')
        gdata.prom(ctr,2) = output.data.sweep{1}.findpeaks.sweepprom(1,1);
        gdata.adaptation2(ctr,:) = output.data.sweep{1}.findpeaks.sweepprom(:,1);
        gdata.pksinitiationvoltage(ctr,2) = output.Vmon.raw(1,output.data.sweep{1}.plots.pksinitiationindex(1,1));
        gdata.timing(ctr,2) = output.data.sweep{1}.plots.locAfterStim(1,1);
        gdata.perdiff(ctr,2) = output.data.sweep{1}.findpeaks.sweepperdiff(1,1);
        gdata.events(ctr,2) = output.data.sweep{1}.findpeaks.sweepevents(1,1);
        gdata.halfwidth(ctr,2) = output.data.sweep{1}.plots.halfwidth(1,1);
        gdata.fnames2{ctr,1} = filesgroup2(ctr).name;
        gdata.Imon.raw2{ctr} = output.Imon.raw(1,:);
        gdata.time2{ctr} = output.data.sweep{1}.time{1};
        gdata.integral(ctr,2) = output.data.sweep{1}.plots.integral(1,1);
        gdata.onsettime(ctr,2) = output.data.sweep{1}.plots.onsettime(1,1);
        gdata.offsettime(ctr,2) = output.data.sweep{1}.plots.offsettime(1,1);
        gdata.onsetcurrent(ctr,2) = output.data.sweep{1}.plots.onsetcurrent(1,1);
        gdata.offsetcurrent(ctr,2) = output.data.sweep{1}.plots.offsetcurrent(1,1);
        gdata.onsetslopecurrent(ctr,2) = output.data.sweep{1}.plots.onsetslopecurrent(1,1);
        gdata.offsetslopecurrent(ctr,2) = output.data.sweep{1}.plots.offsetslopecurrent(1,1);
    end
    gdata.prom(gdata.prom == 0) = NaN;
    gdata.pksinitiationvoltage(gdata.pksinitiationvoltage == 0) = NaN;
    gdata.timing(gdata.timing == 0) = NaN;
    gdata.events(gdata.events == 0) = NaN;
    gdata.perdiff(gdata.perdiff == 0) = NaN;
    gdata.halfwidth(gdata.halfwidth == 0) = NaN;
    gdata.integral(gdata.integral == 0) = NaN;
    gdata.onsettime(gdata.onsettime == 0) = NaN;
    gdata.offsettime(gdata.offsettime == 0) = NaN;
    gdata.onsetcurrent(gdata.onsetcurrent == 0) = NaN;
    gdata.offsetcurrent(gdata.offsetcurrent == 0) = NaN;
    gdata.onsetslopecurrent(gdata.onsetslopecurrent == 0) = NaN;
    gdata.offsetslopecurrent(gdata.offsetslopecurrent == 0) = NaN;
end

%% plotting

if nargin == 1 % Single group analysis
    
    line_width = 1.5;
    marker_size = 10;
    
    for ctr = 1 : size(filesgroup1,1)
        
        h1 = subplot(4,3,1); hold on;
        plot(gdata.prom(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Peak Prominence','FontSize',12); xlabel('Tooth'); ylabel('Amplitude [nA]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h2 = subplot(4,3,2); hold on;
        plot(gdata.pksinitiationvoltage(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Initiation Voltage','FontSize',12); xlabel('Tooth'); ylabel('Voltage [mV]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h3 = subplot(4,3,3); hold on;
        plot(gdata.timing(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Timing after Stimulus','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h4 = subplot(4,3,4); hold on;
        plot(gdata.perdiff(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('1/2 Peak Prominence Difference','FontSize',12); xlabel('Tooth'); ylabel('Percentage');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h5 = subplot(4,3,5); hold on;
        plot(gdata.events(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('# events','FontSize',12); xlabel('Tooth'); ylabel('#');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h6 = subplot(4,3,6); hold on;
        plot(gdata.halfwidth(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Halfwidth','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h7 = subplot(4,3,7); hold on;
        plot(gdata.integral(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Integral','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h8 = subplot(4,3,8); hold on;
        plot(gdata.onsettime(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('onsettime','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h9 = subplot(4,3,9); hold on;
        plot(gdata.offsettime(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('offsettime','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h10 = subplot(4,3,10); hold on;
        plot(gdata.onsetslopecurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Onsetslopecurrent','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
        
        h11 = subplot(4,3,11); hold on;
        plot(gdata.offsetslopecurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Offsetslopecurrent','FontSize',12); xlabel('Tooth'); ylabel('Time [s]');
        xlim([0 5]); ylim('auto')
        set(gca,'xtick',[1 2 3 4 5])
    end
    
elseif nargin == 4 % Combined group analysis
    
    if strcmp(mode, 'Yes')
        [gdata.hypo(1),gdata.p(1)] = ttest(gdata.prom(:,1),gdata.prom(:,2));
        [gdata.hypo(2),gdata.p(2)] = ttest(gdata.pksinitiationvoltage(:,1),gdata.pksinitiationvoltage(:,2));
        [gdata.hypo(3),gdata.p(3)] = ttest(gdata.timing(:,1),gdata.timing(:,2));
        [gdata.hypo(4),gdata.p(4)] = ttest(gdata.perdiff(:,1),gdata.perdiff(:,2));
        [gdata.hypo(5),gdata.p(5)] = ttest(gdata.halfwidth(:,1),gdata.halfwidth(:,2));
        plotstyle = '-ob';
    else
        [gdata.hypo(1),gdata.p(1)] = ttest2(gdata.prom(:,1),gdata.prom(:,2));
        [gdata.hypo(2),gdata.p(2)] = ttest2(gdata.pksinitiationvoltage(:,1),gdata.pksinitiationvoltage(:,2));
        [gdata.hypo(3),gdata.p(3)] = ttest2(gdata.timing(:,1),gdata.timing(:,2));
        [gdata.hypo(4),gdata.p(4)] = ttest2(gdata.perdiff(:,1),gdata.perdiff(:,2));
        [gdata.hypo(5),gdata.p(5)] = ttest2(gdata.halfwidth(:,1),gdata.halfwidth(:,2));
        plotstyle = 'ob';
    end
    
    line_width = 1.5;
    marker_size = 10;
    condition = {group1,group2};
    
    figure();
    set(gcf,'position', get(0,'screensize'));
    
    gdata.prom(gdata.prom == 0) = NaN;
    gdata.pksinitiationvoltage(gdata.pksinitiationvoltage == 0) = NaN;
    gdata.timing(gdata.timing == 0) = NaN;
    gdata.events(gdata.events == 0) = NaN;
    gdata.perdiff(gdata.perdiff == 0) = NaN;
    gdata.halfwidth(gdata.halfwidth == 0) = NaN;
    
    % Fig1
    for ctr = 1 : size(gdata.prom,1)
        h1 = subplot(3,2,1); hold on;
        plot(gdata.prom(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Amplitude. p = ' num2str(gdata.p(1))],'FontSize',16); xlabel('Condition'); ylabel('Amplitude [nA]','FontSize',16);
        
        h2 = subplot(3,2,2); hold on;
        plot(gdata.pksinitiationvoltage(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['Initiation Voltage. p = ' num2str(gdata.p(2))],'FontSize',16); xlabel('Condition'); ylabel('Voltage [mV]','FontSize',16);
        
        h3 = subplot(3,2,3); hold on;
        plot(gdata.timing(ctr,:).*1000, plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Latency. p = ' num2str(gdata.p(3))],'FontSize',16); xlabel('Condition'); ylabel('Time [ms]','FontSize',16);
        
        h4 = subplot(3,2,4); hold on;
        plot(gdata.perdiff(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1-2 Peak Prominence Difference. p = ' num2str(gdata.p(4))],'FontSize',16); xlabel('Condition'); ylabel('Percentage','FontSize',16);
        
        h5 = subplot(3,2,5); hold on;
        plot(gdata.events(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('# events','FontSize',16); xlabel('Condition','FontSize',16); ylabel('#');
        
        h6 = subplot(3,2,6); hold on;
        plot(gdata.halfwidth(ctr,:).*1000, plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title(['1st Peak Halfwidth. p = ' num2str(gdata.p(5))],'FontSize',16); xlabel('Condition'); ylabel('Time [ms]','FontSize',16);
    end
    
    boxplot_position = [ones(size(gdata.prom,1), 1)*0.5, ones(size(gdata.prom,1), 1)*2.5];
    boxplot_width = [0.3, 0.3];
    
    axes(h1); hold on
    boxplot(gdata.prom(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.prom(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h2); hold on
    boxplot(gdata.pksinitiationvoltage(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.pksinitiationvoltage(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h3); hold on
    boxplot(gdata.timing(:,1).*1000 , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.timing(:,2).*1000 , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h4); hold on
    boxplot(gdata.perdiff(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.perdiff(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h5); hold on
    boxplot(gdata.events(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.events(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h6); hold on
    boxplot(gdata.halfwidth(:,1).*1000 , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.halfwidth(:,2).*1000 , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    save(['Groupanalysis ' filename], 'gdata')
    set(gcf,'Position', get(0, 'Screensize'));
    saveas(gcf,[filename ' Boxplots1.png']);
    savefig([filename ' Boxplots1']);
    close all
    
    % Fig 2
    for ctr = 1 : size(gdata.prom,1)
        h7 = subplot(3,2,1); hold on;
        plot(gdata.onsettime(ctr,:).*1000, 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Onset Duration','FontSize',12); xlabel('Condition','FontSize',16); ylabel('ms');
        
        h8 = subplot(3,2,2); hold on;
        plot(gdata.offsettime(ctr,:).*1000, 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Offset Duration','FontSize',12); xlabel('Condition','FontSize',16); ylabel('ms');
        
        h9 = subplot(3,2,3); hold on;
        plot(gdata.onsetcurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Onset Current','FontSize',12); xlabel('Condition','FontSize',16); ylabel('nA');
        
        h10 = subplot(3,2,4); hold on;
        plot(gdata.offsetcurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Offset Current','FontSize',12); xlabel('Condition','FontSize',16); ylabel('nA');
        
        h11 = subplot(3,2,5); hold on;
        plot(gdata.onsetslopecurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Onset Current/Duration','FontSize',12); xlabel('Condition','FontSize',16); ylabel('nA/ms');
        
        h12 = subplot(3,2,6); hold on;
        plot(gdata.offsetslopecurrent(ctr,:), 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Offset Current/Duration','FontSize',12); xlabel('Condition','FontSize',16); ylabel('nA/ms');
    end
    
    axes(h7); hold on
    boxplot(gdata.onsettime(:,1).*1000 , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.onsettime(:,2).*1000 , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h8); hold on
    boxplot(gdata.offsettime(:,1).*1000 , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.offsettime(:,2).*1000 , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h9); hold on
    boxplot(gdata.onsetcurrent(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.onsetcurrent(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h10); hold on
    boxplot(gdata.offsetcurrent(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.offsetcurrent(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h11); hold on
    boxplot(gdata.onsetslopecurrent(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.onsetslopecurrent(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    axes(h12); hold on
    boxplot(gdata.offsetslopecurrent(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.offsetslopecurrent(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    save(['Groupanalysis ' filename], 'gdata')
    set(gcf,'Position', get(0, 'Screensize'));
    saveas(gcf,[filename ' Boxplots2.png']);
    savefig([filename ' Boxplots2']);
    close all
    
    % fig 3
    for ctr = 1 : size(gdata.prom,1)
        h13 = subplot(3,2,1); hold on;
        plot(gdata.integral(ctr,:).*1000, 'ob', 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
        title('Integral','FontSize',12); xlabel('Condition','FontSize',16); ylabel('Time [ms]');
    end
    
    axes(h13); hold on
    boxplot(gdata.integral(:,1).*1000 , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
    boxplot(gdata.integral(:,2).*1000 , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
    xlim([0 3]); ylim('auto')
    set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
    
    save(['Groupanalysis ' filename], 'gdata')
    set(gcf,'Position', get(0, 'Screensize'));
    saveas(gcf,[filename ' Boxplots3.png']);
    savefig([filename ' Boxplots3']);
    close all
end
%% finding the most representable experiments.
if nargin == 4
    ValuesNormalizedProm1=(gdata.prom(:,1))./(max(gdata.prom(:,1)));
    ValuesNormalizedInivolt1=(gdata.pksinitiationvoltage(:,1))./(min(gdata.pksinitiationvoltage(:,1)));
    ValuesNormalizedTiming1=(gdata.timing(:,1))./(max(gdata.timing(:,1)));
    ValuesNormalizedPeakdif1=(gdata.perdiff(:,1))./(max(gdata.perdiff(:,1)));
    ValuesNormalizedHalfwidth1=(gdata.halfwidth(:,1))./(max(gdata.halfwidth(:,1)));
    ValuesNormalizedEvent1=(gdata.events(:,1))./(max(gdata.events(:,1)));
    
    ValuesNormalizedProm2=(gdata.prom(:,2))./(max(gdata.prom(:,2)));
    ValuesNormalizedInivolt2=(gdata.pksinitiationvoltage(:,2))./(min(gdata.pksinitiationvoltage(:,2)));
    ValuesNormalizedTiming2=(gdata.timing(:,2))./(max(gdata.timing(:,2)));
    ValuesNormalizedPeakdif2=(gdata.perdiff(:,2))./(max(gdata.perdiff(:,2)));
    ValuesNormalizedHalfwidth2=(gdata.halfwidth(:,2))./(max(gdata.halfwidth(:,2)));
    ValuesNormalizedEvent2=(gdata.events(:,2))./(max(gdata.events(:,2)));
    
    NormalizedMeanProm1 = nanmean(ValuesNormalizedProm1,1);
    NormalizedMeanInivolt1 = nanmean(ValuesNormalizedInivolt1,1);
    NormalizedMeanTiming1 = nanmean(ValuesNormalizedTiming1,1);
    NormalizedMeanPeakdif1 = nanmean(ValuesNormalizedPeakdif1,1);
    NormalizedMeanHalfwidth1 = nanmean(ValuesNormalizedHalfwidth1,1);
    NormalizedMeanEvent1 = nanmean(ValuesNormalizedEvent1,1);
    
    NormalizedMeanProm2 = nanmean(ValuesNormalizedProm2,1);
    NormalizedMeanInivolt2 = nanmean(ValuesNormalizedInivolt2,1);
    NormalizedMeanTiming2 = nanmean(ValuesNormalizedTiming2,1);
    NormalizedMeanPeakdif2 = nanmean(ValuesNormalizedPeakdif2,1);
    NormalizedMeanHalfwidth2 = nanmean(ValuesNormalizedHalfwidth2,1);
    NormalizedMeanEvent2 = nanmean(ValuesNormalizedEvent2,1);
    
    difference1 = zeros(size(gdata.fnames1,1),6);
    difference2 = zeros(size(gdata.fnames2,1),6);
    
    for i=1:size(gdata.fnames1,1)
        difference1(i,1) = abs(NormalizedMeanProm1-ValuesNormalizedProm1(i,:));
        difference1(i,2) = abs(NormalizedMeanInivolt1-ValuesNormalizedInivolt1(i,:));
        difference1(i,3) = abs(NormalizedMeanTiming1-ValuesNormalizedTiming1(i,:));
        difference1(i,4) = abs(NormalizedMeanPeakdif1-ValuesNormalizedPeakdif1(i,:));
        difference1(i,5) = abs(NormalizedMeanHalfwidth1-ValuesNormalizedHalfwidth1(i,:));
        difference1(i,6) = abs(NormalizedMeanEvent1-ValuesNormalizedEvent1(i,:));
    end
    for i=1:size(gdata.fnames2,1)
        difference2(i,1) = abs(NormalizedMeanProm2-ValuesNormalizedProm2(i,:));
        difference2(i,2) = abs(NormalizedMeanInivolt2-ValuesNormalizedInivolt2(i,:));
        difference2(i,3) = abs(NormalizedMeanTiming2-ValuesNormalizedTiming2(i,:));
        difference2(i,4) = abs(NormalizedMeanPeakdif2-ValuesNormalizedPeakdif2(i,:));
        difference2(i,5) = abs(NormalizedMeanHalfwidth2-ValuesNormalizedHalfwidth2(i,:));
        difference2(i,6) = abs(NormalizedMeanEvent2-ValuesNormalizedEvent2(i,:));
    end
    
    CumulativeDifference1 = sum(difference1,2);
    CumulativeDifference2 = sum(difference2,2);
    
    optimal1 = find(CumulativeDifference1 == min(CumulativeDifference1));
    optimal2 = find(CumulativeDifference2 == min(CumulativeDifference2));
    
    figure(2)
    subplot(2,1,1)
    plot(gdata.time1{optimal1},gdata.Imon.raw1{optimal1},'b')
    title([group1 ' representation. N = ' num2str(size(filesgroup1,1))],'FontSize',16)
    ylabel('Current (nA)','FontSize',16)
    xlabel('Time (s)','FontSize',16)
    
    subplot(2,1,2)
    plot(gdata.time2{optimal2},gdata.Imon.raw2{optimal2},'b')
    title([group2 ' representation. N = ' num2str(size(filesgroup2,1))],'FontSize',16)
    ylabel('Current (nA)','FontSize',16)
    xlabel('Time (s)','FontSize',16)
end

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,[filename ' Representation plots.png']);
savefig([filename ' Representation plots']);
close all

%% Normalized plots
% 
% if strcmp(mode, 'Yes')
%     gdata.promnorm = gdata.prom./gdata.prom(:,1);
%     gdata.pksinitiationvoltagenorm = gdata.pksinitiationvoltage./gdata.pksinitiationvoltage(:,1);
%     gdata.timingnorm = gdata.timing./gdata.timing(:,1);
%     gdata.perdiffnorm = gdata.perdiff./gdata.perdiff(:,1);
%     gdata.halfwidthnorm = gdata.halfwidth./gdata.halfwidth(:,1);
%     
%     [gdata.hypon(1),gdata.pn(1)] = ttest2(gdata.promnorm(:,1),gdata.promnorm(:,2));
%     [gdata.hypon(2),gdata.pn(2)] = ttest2(gdata.pksinitiationvoltagenorm(:,1),gdata.pksinitiationvoltagenorm(:,2));
%     [gdata.hypon(3),gdata.pn(3)] = ttest2(gdata.timingnorm(:,1),gdata.timingnorm(:,2));
%     [gdata.hypon(4),gdata.pn(4)] = ttest2(gdata.perdiffnorm(:,1),gdata.perdiffnorm(:,2));
%     [gdata.hypon(5),gdata.pn(5)] = ttest2(gdata.halfwidthnorm(:,1),gdata.halfwidthnorm(:,2));
%     
%     line_width = 1.5;
%     marker_size = 10;
%     condition = {group1,group2};
%     
%     figure();
%     set(gcf,'position', get(0,'screensize'));
%     
%     gdata.prom(gdata.prom == 0) = NaN;
%     gdata.pksinitiationvoltage(gdata.pksinitiationvoltage == 0) = NaN;
%     gdata.timing(gdata.timing == 0) = NaN;
%     gdata.events(gdata.events == 0) = NaN;
%     gdata.perdiff(gdata.perdiff == 0) = NaN;
%     gdata.halfwidth(gdata.halfwidth == 0) = NaN;
%     
%     for ctr = 1 : size(gdata.prom,1)
%         h1 = subplot(3,2,1); hold on;
%         plot(gdata.promnorm(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title(['1st Peak Amplitude. p = ' num2str(gdata.pn(1))],'FontSize',16); xlabel('Condition'); ylabel('Amplitude [nA]','FontSize',16);
%         
%         h2 = subplot(3,2,2); hold on;
%         plot(gdata.pksinitiationvoltagenorm(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title(['Initiation Voltage. p = ' num2str(gdata.pn(2))],'FontSize',16); xlabel('Condition'); ylabel('Voltage [mV]','FontSize',16);
%         
%         h3 = subplot(3,2,3); hold on;
%         plot(gdata.timingnorm(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title(['1st Peak Latency. p = ' num2str(gdata.pn(3))],'FontSize',16); xlabel('Condition'); ylabel('Time [s]','FontSize',16);
%         
%         h4 = subplot(3,2,4); hold on;
%         plot(gdata.perdiffnorm(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title(['1-2 Peak Prominence Difference. p = ' num2str(gdata.pn(4))],'FontSize',16); xlabel('Condition'); ylabel('Percentage','FontSize',16);
%         
%         h5 = subplot(3,2,5); hold on;
%         plot(gdata.events(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title('# events','FontSize',16); xlabel('Condition','FontSize',16); ylabel('#');
%         
%         h6 = subplot(3,2,6); hold on;
%         plot(gdata.halfwidthnorm(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
%         title(['1st Peak Halfwidth. p = ' num2str(gdata.pn(5))],'FontSize',16); xlabel('Condition'); ylabel('Time [s]','FontSize',16);
%     end
%     
%     
%     boxplot_position = [ones(size(gdata.prom,1), 1)*0.5, ones(size(gdata.prom,1), 1)*2.5];
%     boxplot_width = [0.3, 0.3];
%     
%     axes(h1); hold on
%     boxplot(gdata.promnorm(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.promnorm(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     axes(h2); hold on
%     boxplot(gdata.pksinitiationvoltagenorm(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.pksinitiationvoltagenorm(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     axes(h3); hold on
%     boxplot(gdata.timingnorm(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.timingnorm(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     axes(h4); hold on
%     boxplot(gdata.perdiffnorm(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.perdiffnorm(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     axes(h5); hold on
%     boxplot(gdata.events(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.events(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     axes(h6); hold on
%     boxplot(gdata.halfwidthnorm(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
%     boxplot(gdata.halfwidthnorm(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
%     xlim([0 3]); ylim('auto')
%     set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)
%     
%     set(gcf,'Position', get(0, 'Screensize'));
%     saveas(gcf,'normalized D1 test Boxplots.png');
%     savefig('normalized D1 test Boxplots');
%     close all
% end











