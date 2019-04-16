load('Groupanalysis E vs I 10ms.mat')
gdata10 = gdata;
load('Groupanalysis E vs I 50ms.mat')
gdata50 = gdata;
load('Groupanalysis E vs I 100ms.mat')
gdata100 = gdata;

clear gdata

cdata.exc.adaptation100(1,1:5) = nanmean(gdata100.excadap,1);
cdata.inh.adaptation100(1,1:5) = nanmean(gdata100.inhadap,1);
cdata.exc.adaptation50(1,1:5) = nanmean(gdata50.excadap,1);
cdata.inh.adaptation50(1,1:5) = nanmean(gdata50.inhadap,1);
cdata.exc.adaptation10(1,1:5) = nanmean(gdata10.excadap,1);
cdata.inh.adaptation10(1,1:5) = nanmean(gdata10.inhadap,1);
% 
% cdata.exc.adaptation100norm = nanmean(gdata100.adaptation1./gdata100.adaptation1(:,1));
% cdata.inh.adaptation100norm = nanmean(gdata100.adaptation2./gdata100.adaptation2(:,1));
% cdata.exc.adaptation50norm = nanmean(gdata50.adaptation1./gdata50.adaptation1(:,1));
% cdata.inh.adaptation50norm = nanmean(gdata50.adaptation2./gdata50.adaptation2(:,1));
% cdata.exc.adaptation10norm = nanmean(gdata10.adaptation1./gdata10.adaptation1(:,1));
% cdata.inh.adaptation10norm = nanmean(gdata10.adaptation2./gdata10.adaptation2(:,1));

cdata.exc.adaptation100std = nanstd(gdata100.excadap,1);
cdata.inh.adaptation100std = nanstd(gdata100.inhadap,1);
cdata.exc.adaptation50std = nanstd(gdata50.excadap,1);
cdata.inh.adaptation50std = nanstd(gdata50.inhadap,1);
cdata.exc.adaptation10std = nanstd(gdata10.excadap,1);
cdata.inh.adaptation10std = nanstd(gdata10.inhadap,1);

gdata10.offset(gdata10.offset == 0) = nan;
gdata10.onset(gdata10.onset == 0) = nan;
gdata10.halfwidth(gdata10.halfwidth == 0) = nan;

%% p values

h(1) = ABKS(gdata100.halfwidth(:,1));
h(2) = ABKS(gdata100.halfwidth(:,2));
h(3) = ABKS(gdata50.halfwidth(:,1));
h(4) = ABKS(gdata50.halfwidth(:,2));
h(5) = ABKS(gdata10.halfwidth(:,1));
h(6) = ABKS(gdata10.halfwidth(:,2));
h(7) = ABKS(gdata100.onset(:,1));
h(8) = ABKS(gdata100.onset(:,2));
h(9) = ABKS(gdata50.onset(:,1));
h(10) = ABKS(gdata50.onset(:,2));
h(11) = ABKS(gdata10.onset(:,1));
h(12) = ABKS(gdata10.onset(:,2));
h(13) = ABKS(gdata100.offset(:,1));
h(14) = ABKS(gdata100.offset(:,2));
h(15) = ABKS(gdata50.offset(:,1));
h(16) = ABKS(gdata50.offset(:,2));
h(17) = ABKS(gdata10.offset(:,1));
h(18) = ABKS(gdata10.offset(:,2));

[~,p(1)] = ttest(gdata100.halfwidth(:,1),gdata100.halfwidth(:,2));
[~,p(2)] = ttest(gdata50.halfwidth(:,1),gdata50.halfwidth(:,2));
[~,p(3)] = ttest(gdata10.halfwidth(:,1),gdata10.halfwidth(:,2));
[~,p(4)] = ttest(gdata100.onset(:,1),gdata100.onset(:,2));
[~,p(5)] = ttest(gdata50.onset(:,1),gdata50.onset(:,2));
[~,p(6)] = ttest(gdata10.onset(:,1),gdata10.onset(:,2));
[~,p(7)] = ttest(gdata100.offset(:,1),gdata100.offset(:,2));
[~,p(8)] = ttest(gdata50.offset(:,1),gdata50.offset(:,2));
[~,p(9)] = ttest(gdata10.offset(:,1),gdata10.offset(:,2));

%%
subplot(4,3,1)
hold on ;grid on
plot(1,gdata100.halfwidth(:,1)*1000,'ro')
plot(2,gdata100.halfwidth(:,2)*1000,'bo')
boxplot(gdata100.halfwidth(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata100.halfwidth(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
ylabel('Halfwidth (ms)')
title(['P = ' num2str(p(1))])
xlim([0 3])
ylim([0 2.5])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,2)
hold on ;grid on
plot(1,gdata50.halfwidth(:,1)*1000,'ro')
plot(2,gdata50.halfwidth(:,2)*1000,'bo')
boxplot(gdata50.halfwidth(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata50.halfwidth(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(2))])
xlim([0 3])
ylim([0 2.5])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,3)
hold on ;grid on
plot(1,gdata10.halfwidth(:,1)*1000,'ro')
plot(2,gdata10.halfwidth(:,2)*1000,'bo')
boxplot(gdata10.halfwidth(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata10.halfwidth(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(3))])
xlim([0 3])
ylim([0 1.55])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,4)
hold on ;grid on
plot(1,gdata100.onset(:,1)*1000,'ro')
plot(2,gdata100.onset(:,2)*1000,'bo')
boxplot(gdata100.onset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata100.onset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(4))])
ylabel('Onsettime (ms)')
xlim([0 3])
ylim([0 4])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,5)
hold on ;grid on
plot(1,gdata50.onset(:,1)*1000,'ro')
plot(2,gdata50.onset(:,2)*1000,'bo')
boxplot(gdata50.onset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata50.onset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(5))])
xlim([0 3])
ylim([0 4])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,6)
hold on ;grid on
plot(1,gdata10.onset(:,1)*1000,'ro')
plot(2,gdata10.onset(:,2)*1000,'bo')
boxplot(gdata10.onset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata10.onset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(6))])
xlim([0 3])
ylim([0 2])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,7)
hold on ;grid on
plot(1,gdata100.offset(:,1)*1000,'ro')
plot(2,gdata100.offset(:,2)*1000,'bo')
boxplot(gdata100.offset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata100.offset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
ylabel('Offsettime (ms)')
title(['P = ' num2str(p(7))])
xlim([0 3])
ylim([0 6])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,8)
hold on ;grid on
plot(1,gdata50.offset(:,1)*1000,'ro')
plot(2,gdata50.offset(:,2)*1000,'bo')
boxplot(gdata50.offset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata50.offset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(8))])
xlim([0 3])
ylim([0 6])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,9)
hold on ;grid on
plot(1,gdata10.offset(:,1)*1000,'ro')
plot(2,gdata10.offset(:,2)*1000,'bo')
boxplot(gdata10.offset(:,1)*1000, 'MedianStyle', 'line', 'Positions', 0.5);
boxplot(gdata10.offset(:,2)*1000, 'MedianStyle', 'line', 'Positions', 2.5);
title(['P = ' num2str(p(9))])
xlim([0 3])
ylim([0 6])
set(gca,'xtick',[1 2],'xticklabel',{'Excitatory','Inhibitory'},'FontSize',12)

subplot(4,3,10)
hold on; grid on
pl1 = errorbar(0.9:1:4.9,cdata.exc.adaptation100,cdata.exc.adaptation100std,'ro-');
pl2 = errorbar(1.1:1:5.1,cdata.inh.adaptation100,cdata.inh.adaptation100std,'bo-');
ylim([0 1.2])
xlabel('Sawtooth nr','FontSize',12)
xticks([1 2 3 4 5])
ylabel('Normalized adaptation amplitude','FontSize',12)

subplot(4,3,11)
hold on; grid on
pl1 = errorbar(0.9:1:4.9,cdata.exc.adaptation50,cdata.exc.adaptation50std,'ro-');
pl2 = errorbar(1.1:1:5.1,cdata.inh.adaptation50,cdata.inh.adaptation50std,'bo-');
ylim([0 1.2])
xlabel('Sawtooth nr','FontSize',16)
xticks([1 2 3 4 5])

subplot(4,3,12)
hold on; grid on
pl1 = errorbar(0.9:1:4.9,cdata.exc.adaptation10,cdata.exc.adaptation10std,'ro-');
pl2 = errorbar(1.1:1:5.1,cdata.inh.adaptation10,cdata.inh.adaptation10std,'bo-');
ylim([0 1.2])
xlabel('Sawtooth nr','FontSize',16)
xticks([1 2 3 4 5])

%% Within sweep amplitude adaptation

figure(1); 
set(gcf,'position', get(0,'screensize'));
hold on; grid on; box on

x = 1:5;

subplot(1,3,1); hold on; grid on; box on
p1 = errorbar(cdata.exc.adaptation100norm,cdata.exc.adaptation100std,'ro-');
p2 = errorbar(cdata.inh.adaptation100norm,cdata.inh.adaptation100std,'bo-');
ylim([0 1.5])
title('100ms (10 Hz)','FontSize',22)
xlabel('Tooth','FontSize',16)
xticks([1 2 3 4 5])
ylabel('Normalized amplitude','FontSize',16)
legend([p1 p2],{'Excitatory','Inhibitory'})

subplot(1,3,2); hold on; grid on; box on
errorbar(cdata.exc.adaptation50norm,cdata.exc.adaptation50std,'ro-');
errorbar(cdata.inh.adaptation50norm,cdata.inh.adaptation50std,'bo-');
ylim([0 1.5])
title('50ms (20 Hz)','FontSize',22)
xlabel('Tooth','FontSize',16)
xticks([1 2 3 4 5])
ylabel('Normalized amplitude','FontSize',16)

subplot(1,3,3); hold on; grid on; box on
errorbar(cdata.exc.adaptation10norm,cdata.exc.adaptation10std,'ro-');
errorbar(cdata.inh.adaptation10norm,cdata.inh.adaptation10std,'bo-');
ylim([0 1.5])
title('10ms (100 Hz)','FontSize',22)
xlabel('Tooth','FontSize',16)
xticks([1 2 3 4 5])
ylabel('Normalized amplitude','FontSize',16)

saveas(gcf,'VCsaw_Within sweep adaptation.png');
savefig('VCsaw_Within sweep adaptation');

%% Combined parameters

line_width = 1.5;
marker_size = 10;
boxplot_position1 = [ones(size(gdata100.prom,1), 1)*0.5, ones(size(gdata100.prom,1), 1)*2.5];
boxplot_position2 = [ones(size(gdata50.prom,1), 1)*3.5, ones(size(gdata50.prom,1), 1)*5.5];
boxplot_position3 = [ones(size(gdata10.prom,1), 1)*6.5, ones(size(gdata10.prom,1), 1)*8.5];
boxplot_position4 = [ones(size(gdata100.prom,1), 1)*1, ones(size(gdata100.prom,1), 1)*4.5];
boxplot_position5 = [ones(size(gdata50.prom,1), 1)*2, ones(size(gdata50.prom,1), 1)*5.5];
boxplot_position6 = [ones(size(gdata10.prom,1), 1)*3, ones(size(gdata10.prom,1), 1)*6.5];
boxplot_width = [0.3, 0.3];
plotstyle = {'ob','or','og'};
condition = {'Excitatory','Inhibitory','Excitatory','Inhibitory','Excitatory','Inhibitory'};

gdata100.prom(gdata100.prom == 0) = NaN;
gdata100.pksinitiationvoltage(gdata100.pksinitiationvoltage == 0) = NaN;
gdata100.timing(gdata100.timing == 0) = NaN;
gdata100.events(gdata100.events == 0) = NaN;
gdata100.perdiff(gdata100.perdiff == 0) = NaN;
gdata100.halfwidth(gdata100.halfwidth == 0) = NaN;

%%
figure(2); 
set(gcf,'position', get(0,'screensize'));
hold on

for ctr = 1 : size(gdata100.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([1 2],gdata100.prom(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
  
    h2 = subplot(3,1,2); hold on;
    plot([1 2],gdata100.pksinitiationvoltage(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
   
    h3 = subplot(3,1,3); hold on;
    plot([1 2],gdata100.timing(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata50.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([4 5],gdata50.prom(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
   
    h2 = subplot(3,1,2); hold on;
    plot([4 5],gdata50.pksinitiationvoltage(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h3 = subplot(3,1,3); hold on;
    plot([4 5],gdata50.timing(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata10.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([7 8],gdata10.prom(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Amplitude','FontSize',16); ylabel('Amplitude [nA]','FontSize',16);
    
    h2 = subplot(3,1,2); hold on;
    plot([7 8],gdata10.pksinitiationvoltage(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('Initiation Voltage','FontSize',16); ylabel('Voltage [mV]','FontSize',16);
    
    h3 = subplot(3,1,3); hold on;
    plot([7 8],gdata10.timing(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Latency','FontSize',16); ylabel('Time [s]','FontSize',16);
end

axes(h1); hold on
boxplot(gdata100.prom(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.prom(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.prom(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.prom(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.prom(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.prom(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

axes(h2); hold on
boxplot(gdata100.pksinitiationvoltage(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.pksinitiationvoltage(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.pksinitiationvoltage(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.pksinitiationvoltage(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.pksinitiationvoltage(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.pksinitiationvoltage(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

axes(h3); hold on
boxplot(gdata100.timing(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.timing(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.timing(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.timing(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.timing(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.timing(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

suptitle('Sawtooth parameter comparison: 10 Hz (Blue), 20 Hz (red), 100 Hz (green)')

saveas(gcf,'VCsaw_I vs E + frequencies parameters (Part 1 Version 1).png');
savefig('VCsaw_I vs E + frequencies parameters (Part 1 Version 1)');

%%
figure(3)
set(gcf,'position', get(0,'screensize'));
hold on

for ctr = 1 : size(gdata100.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([1 2],gdata100.perdiff(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h5 = subplot(3,1,2); hold on;
    plot([1 2],gdata100.events(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h6 = subplot(3,1,3); hold on;
    plot([1 2],gdata100.halfwidth(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata50.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([4 5],gdata50.perdiff(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h5 = subplot(3,1,2); hold on;
    plot([4 5],gdata50.events(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h6 = subplot(3,1,3); hold on;
    plot([4 5],gdata50.halfwidth(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata10.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([7 8],gdata10.perdiff(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1-2 Peak Prominence Difference','FontSize',16); ylabel('Percentage','FontSize',16);
    
    h5 = subplot(3,1,2); hold on;
    plot([7 8],gdata10.events(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('# events','FontSize',16); ylabel('#');
    
    h6 = subplot(3,1,3); hold on;
    plot([7 8],gdata10.halfwidth(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Halfwidth','FontSize',16); ylabel('Time [s]','FontSize',16);
end

axes(h4); hold on
boxplot(gdata100.perdiff(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.perdiff(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.perdiff(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.perdiff(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.perdiff(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.perdiff(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

axes(h5); hold on
boxplot(gdata100.events(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.events(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.events(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.events(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.events(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.events(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

axes(h6); hold on
boxplot(gdata100.halfwidth(:,1) , boxplot_position1(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata100.halfwidth(:,2) , boxplot_position1(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
boxplot(gdata50.halfwidth(:,1) , boxplot_position2(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3.5 , 'Widths', boxplot_width);
boxplot(gdata50.halfwidth(:,2) , boxplot_position2(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.halfwidth(:,1) , boxplot_position3(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
boxplot(gdata10.halfwidth(:,2) , boxplot_position3(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 8.5 , 'Widths', boxplot_width);
xlim([0 9]); ylim('auto')
set(gca,'xtick',[1 2 4 5 7 8],'xticklabel',condition,'FontSize',16)

suptitle('Sawtooth parameter comparison: 10 Hz (Blue), 20 Hz (red), 100 Hz (green)')

saveas(gcf,'VCsaw_I vs E + frequencies parameters (Part 2 Version 1).png');
savefig('VCsaw_I vs E + frequencies parameters (Part 2 Version 1)');

%%
figure(4)
set(gcf,'position', get(0,'screensize'));
hold on

for ctr = 1 : size(gdata100.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([0.5 4],gdata100.prom(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
  
    h2 = subplot(3,1,2); hold on;
    plot([0.5 4],gdata100.pksinitiationvoltage(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
   
    h3 = subplot(3,1,3); hold on;
    plot([0.5 4],gdata100.timing(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata50.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([1.5 5],gdata50.prom(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
   
    h2 = subplot(3,1,2); hold on;
    plot([1.5 5],gdata50.pksinitiationvoltage(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h3 = subplot(3,1,3); hold on;
    plot([1.5 5],gdata50.timing(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata10.prom,1)
    h1 = subplot(3,1,1); hold on;
    plot([2.5 6],gdata10.prom(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Amplitude','FontSize',16); ylabel('Amplitude [nA]','FontSize',16);
    
    h2 = subplot(3,1,2); hold on;
    plot([2.5 6],gdata10.pksinitiationvoltage(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('Initiation Voltage','FontSize',16); ylabel('Voltage [mV]','FontSize',16);
    
    h3 = subplot(3,1,3); hold on;
    plot([2.5 6],gdata10.timing(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Latency','FontSize',16); ylabel('Time [s]','FontSize',16);
end

axes(h1); hold on
boxplot(gdata100.prom(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1, 'Widths', boxplot_width);
boxplot(gdata100.prom(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5, 'Widths', boxplot_width);
boxplot(gdata50.prom(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2, 'Widths', boxplot_width);
boxplot(gdata50.prom(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5, 'Widths', boxplot_width);
boxplot(gdata10.prom(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3, 'Widths', boxplot_width);
boxplot(gdata10.prom(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5, 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

axes(h2); hold on
boxplot(gdata100.pksinitiationvoltage(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1, 'Widths', boxplot_width);
boxplot(gdata100.pksinitiationvoltage(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5, 'Widths', boxplot_width);
boxplot(gdata50.pksinitiationvoltage(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2, 'Widths', boxplot_width);
boxplot(gdata50.pksinitiationvoltage(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5, 'Widths', boxplot_width);
boxplot(gdata10.pksinitiationvoltage(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3, 'Widths', boxplot_width);
boxplot(gdata10.pksinitiationvoltage(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

axes(h3); hold on
boxplot(gdata100.timing(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1, 'Widths', boxplot_width);
boxplot(gdata100.timing(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5, 'Widths', boxplot_width);
boxplot(gdata50.timing(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2, 'Widths', boxplot_width);
boxplot(gdata50.timing(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5, 'Widths', boxplot_width);
boxplot(gdata10.timing(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3, 'Widths', boxplot_width);
boxplot(gdata10.timing(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5, 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

suptitle('Sawtooth parameter comparison: 10 Hz (Blue), 20 Hz (red), 100 Hz (green)')

saveas(gcf,'VCsaw_I vs E + frequencies parameters (Part 1 Version 2).png');
savefig('VCsaw_I vs E + frequencies parameters (Part 1 Version 2)');

%% 
figure(5)
set(gcf,'position', get(0,'screensize'));
hold on

for ctr = 1 : size(gdata100.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([0.5 4],gdata100.perdiff(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h5 = subplot(3,1,2); hold on;
    plot([0.5 4],gdata100.events(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h6 = subplot(3,1,3); hold on;
    plot([0.5 4],gdata100.halfwidth(ctr,:), plotstyle{1}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata50.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([1.5 5],gdata50.perdiff(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h5 = subplot(3,1,2); hold on;
    plot([1.5 5],gdata50.events(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
    
    h6 = subplot(3,1,3); hold on;
    plot([1.5 5],gdata50.halfwidth(ctr,:), plotstyle{2}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);
end
for ctr = 1 : size(gdata10.prom,1)
    h4 = subplot(3,1,1); hold on;
    plot([2.5 6],gdata10.perdiff(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1-2 Peak Prominence Difference','FontSize',16); ylabel('Percentage','FontSize',16);
    
    h5 = subplot(3,1,2); hold on;
    plot([2.5 6],gdata10.events(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('# events','FontSize',16); ylabel('#');
    
    h6 = subplot(3,1,3); hold on;
    plot([2.5 6],gdata10.halfwidth(ctr,:), plotstyle{3}, 'MarkerSize', marker_size, 'LineWidth', line_width, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.5,0.5,0.5]);
    title('1st Peak Halfwidth','FontSize',16); ylabel('Time [s]','FontSize',16);
end

axes(h4); hold on
boxplot(gdata100.perdiff(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata100.perdiff(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5 , 'Widths', boxplot_width);
boxplot(gdata50.perdiff(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata50.perdiff(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.perdiff(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
boxplot(gdata10.perdiff(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

axes(h5); hold on
boxplot(gdata100.events(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata100.events(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5 , 'Widths', boxplot_width);
boxplot(gdata50.events(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata50.events(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.events(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
boxplot(gdata10.events(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

axes(h6); hold on
boxplot(gdata100.halfwidth(:,1) , boxplot_position4(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata100.halfwidth(:,2) , boxplot_position4(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 4.5 , 'Widths', boxplot_width);
boxplot(gdata50.halfwidth(:,1) , boxplot_position5(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata50.halfwidth(:,2) , boxplot_position5(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 5.5 , 'Widths', boxplot_width);
boxplot(gdata10.halfwidth(:,1) , boxplot_position6(:,1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
boxplot(gdata10.halfwidth(:,2) , boxplot_position6(:,2) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 6.5 , 'Widths', boxplot_width);
xlim([0 7]); ylim('auto')
set(gca,'xtick',[1.75 5.25],'xticklabel',condition,'FontSize',16)

suptitle('Sawtooth parameter comparison: 10 Hz (Blue), 20 Hz (red), 100 Hz (green)')

saveas(gcf,'VCsaw_I vs E + frequencies parameters (Part 2 Version 2).png');
savefig('VCsaw_I vs E + frequencies parameters (Part 2 Version 2)');

%%
close all

%% TTX figure






