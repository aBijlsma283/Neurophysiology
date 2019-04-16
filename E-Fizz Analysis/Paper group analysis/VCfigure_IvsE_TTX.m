filesgroup1 = dir (['analyzed_VCsaw*' '_C.mat']);
filesgroup2 = dir (['analyzed_VCsaw*' '_TTX.mat']);

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
end

for ctr = 1:size(filesgroup2,1)
        load (filesgroup2(ctr).name, 'output')
        gdata.prom(ctr,2) = 0;
        gdata.adaptation2(ctr,:) = [0 0 0 0 0];
        gdata.pksinitiationvoltage(ctr,2) = 0;
        gdata.timing(ctr,2) = 0;
        gdata.perdiff(ctr,2) = 0;
        gdata.events(ctr,2) = 0;
        gdata.halfwidth(ctr,2) = 0;
        gdata.fnames2{ctr,1} = filesgroup2(ctr).name;
        gdata.Imon.raw2{ctr} = output.Imon.raw(1,:);
        gdata.time2{ctr} = output.data.sweep{1}.time{1};
end

line_width = 1.5;
marker_size = 10;
condition = {'C','TTX'};

figure();
set(gcf,'position', get(0,'screensize'));
plotstyle = '-ob';

for ctr = 1 : size(gdata.prom,1)
    h1 = subplot(2,4,5); hold on;
    plot(gdata.prom(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Amplitude','FontSize',16); xlabel('Condition'); ylabel('Amplitude [nA]','FontSize',16);
    
    h2 = subplot(2,4,6); hold on;
    plot(gdata.pksinitiationvoltage(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('Initiation Voltage','FontSize',16); xlabel('Condition'); ylabel('Voltage [mV]','FontSize',16);
    
    h3 = subplot(2,4,7); hold on;
    plot(gdata.timing(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Latency','FontSize',16); xlabel('Condition'); ylabel('Time [s]','FontSize',16);
    
    h5 = subplot(2,4,8); hold on;
    plot(gdata.events(ctr,:), plotstyle, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('# events','FontSize',16); xlabel('Condition','FontSize',16); ylabel('#');
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
boxplot(gdata.timing(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata.timing(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim([0 3]); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)

axes(h5); hold on
boxplot(gdata.events(:,1) , boxplot_position(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata.events(:,2) , boxplot_position(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim([0 3]); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',condition,'FontSize',16)

subplot(2,4,1:2)
plot(gdata.time1{1,9},gdata.Imon.raw1{1,9})
hold on
plot(gdata.time2{1,9},gdata.Imon.raw2{1,9},'r')
xlabel('Time [S]','FontSize',16)
ylabel('Current [nA]','FontSize',16)
title('Inhibitory','FontSize',16)
xlim([0 1.1])

subplot(2,4,3:4)
plot(gdata.time1{1,2},gdata.Imon.raw1{1,2})
hold on
plot(gdata.time2{1,2},gdata.Imon.raw2{1,2},'r')
xlabel('Time [S]','FontSize',16)
ylabel('Current [nA]','FontSize',16)
title('Excitatory','FontSize',16)
xlim([0 1.1])

% AL215 Exc
% AB53 Inh