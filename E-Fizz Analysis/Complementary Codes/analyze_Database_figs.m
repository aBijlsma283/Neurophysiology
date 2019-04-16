%% VCsaw

close all;

uchhhh = dir ('analyzed_VCsaw*.mat');

for ctr = 1:size(uchhhh,1)
    load (uchhhh(ctr).name, 'output')
    gdata.prom(ctr,:) = output.data.sweep{1}.findpeaks.sweepprom(:,1);
    gdata.pksinitiationvoltage(ctr,:) = output.Vmon.raw(1,output.data.sweep{1}.plots.pksinitiationindex(1,:));
    gdata.timing(ctr,:) = output.data.sweep{1}.plots.locAfterStim(1,:);
    gdata.halfwidth(ctr,:) = output.data.sweep{1}.plots.halfwidth(1,:);
    gdata.fnames1{ctr,1} = uchhhh(ctr).name;
    gdata.Imon.raw1{ctr} = output.Imon.raw(1,:);
    gdata.time1{ctr} = output.data.sweep{1}.time{1};
end

subplot(4,5,1:5)
plot(gdata.time1{1},gdata.Imon.raw1{1,1},'k','lineWidth',1.25)
hold on; grid on
legend('100 ms')
ylabel('Current (nA)','FontSize',16)
xlabel('time [s]')
title('VC-Saw Protocol','FontSize',22)
set(gca, 'FontName', 'Times New Roman')
xlim([0 1.1])

subplot(4,5,6)
hold on; grid on
plot(gdata.prom(1,:),'ok', 'MarkerFaceColor', 'k')
title('1st peak amplitude','FontSize', 16)
ylabel({'100ms','Current (nA)'},'FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')
ylim([4 6])

subplot(4,5,7)
hold on; grid on
plot(gdata.timing(1,:),'ok', 'MarkerFaceColor', 'k')
title('1st peak latency','FontSize',16)
ylabel('Time (s)','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,8)
hold on; grid on
plot(gdata.halfwidth(1,:),'ok', 'MarkerFaceColor', 'k')
title('peak halfwidth','FontSize',16)
ylabel('time (s)','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,9)
hold on; grid on
plot(gdata.pksinitiationvoltage(1,:),'ok', 'MarkerFaceColor', 'k')
title('1st peak initiation voltage','FontSize',16)
ylabel('Voltage (mV)','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,10)
hold on; grid on
plot(gdata.prom(1,:)./gdata.prom(1,1),'ok', 'MarkerFaceColor', 'k')
title('1st peak adaptation','FontSize',16)
ylabel('Normalized current','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,11)
hold on; grid on
plot(gdata.prom(3,:),'ok', 'MarkerFaceColor', 'k')
ylabel({'50ms','Current (nA)'},'FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')
yticks([5 6])

subplot(4,5,12)
hold on; grid on
plot(gdata.timing(3,:),'ok', 'MarkerFaceColor', 'k')
ylabel('Time (s)','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,13)
hold on; grid on
plot(gdata.halfwidth(3,:),'ok', 'MarkerFaceColor', 'k')
ylabel('time (s)','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,14)
hold on; grid on
plot(gdata.pksinitiationvoltage(3,:),'ok', 'MarkerFaceColor', 'k')
ylabel('Voltage (mV)','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,15)
hold on; grid on
plot(gdata.prom(3,:)./gdata.prom(3,1),'ok', 'MarkerFaceColor', 'k')
ylabel('Normalized current','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,16)
hold on; grid on
plot(gdata.prom(2,:),'ok', 'MarkerFaceColor', 'k')
ylabel({'10ms','Current (nA)'},'FontSize',16)
xlabel('Tooth','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,17)
hold on; grid on
plot(gdata.timing(2,:),'ok', 'MarkerFaceColor', 'k')
ylabel('Time (s)','FontSize',16)
xlabel('Tooth','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,18)
hold on; grid on
plot(gdata.halfwidth(2,:),'ok', 'MarkerFaceColor', 'k')
ylabel('time (s)','FontSize',16)
xlabel('Tooth','FontSize',16)
xlim([0 5])
A = gca;
A.YAxis.Exponent = 0;
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,19)
hold on; grid on
plot(gdata.pksinitiationvoltage(2,:),'ok', 'MarkerFaceColor', 'k')
ylabel('Voltage (mV)','FontSize',16)
xlabel('Tooth','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

subplot(4,5,20)
hold on; grid on
plot(gdata.prom(2,:)./gdata.prom(2,1),'ok', 'MarkerFaceColor', 'k')
ylabel('Normalized current','FontSize',16)
xlabel('Tooth','FontSize',16)
xlim([0 5])
set(gca, 'FontName', 'Times New Roman')

%% VCstep
    
for ctr = 1:max(output.exp.sweep)
    colorgradient = jet(max(output.exp.sweep));
    subplot(2,5,1:5)
    hold on; grid on
    plot(output.data.time,output.Imon.raw(ctr,:),'color',colorgradient(ctr,:))
end

subplot(2,4,1:4)
ylabel('Current [nA]','FontSize',16)
xlabel('Time [s]','FontSize',16)
title('VC-Step Protocol','FontSize',22)
colormap(jet(max(output.exp.sweep)));
c = colorbar;
caxis([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(max(output.exp.sweep))])
c.Label.String = 'Voltage [mV]';
set(gca, 'FontName', 'Times New Roman')

subplot(2,4,5)
hold on; grid on
plot(output.data.plots.stimulusrange,output.data.plots.pks(:,1) + output.data.plots.pksinitiation(:,1),'ok', 'MarkerFaceColor', 'k')
ylabel('Current [nA]','FontSize',16)
xlabel('Membrane potential [mV]','FontSize',16)
title('Amplitude','FontSize',16)
xlim([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(end)])
set(gca, 'FontName', 'Times New Roman')

subplot(2,4,6)
hold on; grid on
for i = 1:size(output.data.plots.pks,1)
    if isnan(output.data.plots.pksindex(i,1)) == 0
        plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.pksindex(i,1)-output.data.plots.stimbeginloc),'ok', 'MarkerFaceColor', 'k')
    end
end
ylabel('Time [s]','FontSize',16)
xlabel('Membrane potential [mV]','FontSize',16)
title('Latency','FontSize',16)
A = gca;
A.YAxis.Exponent = 0;
xlim([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(end)])
set(gca, 'FontName', 'Times New Roman')

subplot(2,4,7)
hold on; grid on
plot(output.data.plots.stimulusrange,output.data.plots.halfwidth,'ok', 'MarkerFaceColor', 'k')
ylabel('Time [s]','FontSize',16)
xlabel('Membrane potential [mV]','FontSize',16)
title('Half width','FontSize',16)
A = gca;
A.YAxis.Exponent = 0;
xlim([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(end)])
set(gca, 'FontName', 'Times New Roman')

subplot(2,4,8)
hold on; grid on
plot(output.data.plots.stimulusrange,-(output.data.plots.pks(:,1) + output.data.plots.pksinitiation(:,1)),'ok', 'MarkerFaceColor', 'k')
ylabel('Current [nA]','FontSize',16)
xlabel('Membrane potential [mV]','FontSize',16)
title('I/V Curve','FontSize',16)
xlim([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(end)])
set(gca, 'FontName', 'Times New Roman')
