function analyze_Plotter(filename,protocol,datafolder)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_plotter is part of the mother function 'analyze_patch' and plots
% the output data of analyze_ST and analyze_VCstep.
% 
% cd(datafolder)
load ([filename '.mat']);

colorplot = [{'bo'},{'r*'},{'g+'}];
color = [{'b'},{'r'},{'g'}]; % [0.2 0.2 0.2] for publications
figure()

if strcmp(protocol,'VCsaw') == 1 || strcmp(protocol,'N2A') == 1
    for ctr = 1:max(output.exp.sweep)
        subplot(5,5,1:10)
        plot(output.data.sweep{ctr}.time{1},output.Imon.raw(ctr,:),'color',color{ctr},'lineWidth',1.25)
        hold on; grid on
        ylabel('Current (nA)','FontSize',12)
        title('VC-Saw Protocol','FontSize',22)
        legend('Sweep 1','Sweep 2','Sweep 3')
        xlim([0 output.data.sweep{1,1}.time{1,1}(end)])
        
        subplot(5,5,11)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.pks + output.data.sweep{ctr}.plots.pksinitiation(1,:),colorplot{ctr})
        title('1st event amplitude','FontSize', 12)
        ylabel('Current (nA)')
        xlim([0 5])
        
        subplot(5,5,12)
        hold on; grid on
        plot(output.data.sweep{ctr}.findpeaks.sweepprom(:,2),colorplot{ctr})
        title('2nd event amplitude','FontSize', 12)
        ylabel('Current (nA)')
        xlim([0 5])
                
        subplot(5,5,13)
        hold on; grid on
        plot(output.data.sweep{ctr}.findpeaks.sweepperdiff(:,1),colorplot{ctr})
        xlim([0 5])
        ylim([0 100])
        ylabel('Amp change (%)')
        title('1-2 amplitude difference %','FontSize', 12)

        subplot(5,5,14)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.pkssaw,colorplot{ctr})
        ylabel('Current (nA)')
        title('Tooth amplitude','FontSize',12)
        xlim([0 5])
                
        subplot(5,5,15)
        hold on; grid on
        for peak = 1:5
            plot(output.data.sweep{ctr}.findpeaks.sweepevents(1,:),colorplot{ctr})
        end
        title('Event count','FontSize', 12)
        xlim([0 5])
                
        subplot(5,5,16)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.locAfterStim(1,:),colorplot{ctr})
        title('1st event latency','FontSize',12)
        ylabel('Time (s)')
        xlim([0 5])
        A = gca;
        A.YAxis.Exponent = 0;

        subplot(5,5,17)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.locAfterStim2(1,:),colorplot{ctr})
        title('2nd event latency','FontSize',12)
        ylabel('Time (s)')
        xlim([0 5])
        A = gca;
        A.YAxis.Exponent = 0;
        
        subplot(5,5,18)
        hold on; grid on
        plot(output.data.sweep{ctr}.findpeaks.sweeptimediff(1,:),colorplot{ctr})
        title('1-2 latency difference','FontSize', 12)
        ylabel('Time (s)')
        xlim([0 5])
        A = gca;
        A.YAxis.Exponent = 0;
        
        subplot(5,5,19)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.halfwidth',colorplot{ctr})
        title('1st event peak halfwidth','FontSize',12)
        ylabel('time (s)')
        xlim([0 5])
        A = gca;
        A.YAxis.Exponent = 0;
        
        subplot(5,5,20)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.initiation(1,:),colorplot{ctr})
        title('1st event initiation voltage','FontSize',12)
        ylabel('Voltage (mV)')
        xlim([0 5])
        
        subplot(5,5,21)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.integral(1,:),colorplot{ctr})
        title('1st event integral','FontSize',12)
        ylabel('Integral (nA*s)')
        xlabel('Repetition')
        xlim([0 5])
        
        subplot(5,5,22)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.onsettime(1,:),colorplot{ctr})
        title('1st event onsettime','FontSize',12)
        ylabel('Time (s)')
        xlabel('Repetition')
        xlim([0 5])
        
        subplot(5,5,23)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.offsettime(1,:),colorplot{ctr})
        title('1st event offsettime','FontSize',12)
        ylabel('Time (s)')
        xlabel('Repetition')
        xlim([0 5])
        
        subplot(5,5,24)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.onsetslopecurrent(1,:),colorplot{ctr})
        title('Onset Current/Duration','FontSize',12)
        ylabel('nA/ms')
        xlabel('Repetition')
        xlim([0 5])
        
        subplot(5,5,25)
        hold on; grid on
        plot(output.data.sweep{ctr}.plots.offsetslopecurrent(1,:),colorplot{ctr})
        title('Offset Current/Duration','FontSize',12)
        ylabel('nA/ms')
        xlabel('Repetition')
        xlim([0 5])
    end

elseif strcmp(protocol,'VCstep') == 1
    
    for ctr = 1:max(output.exp.sweep)
        colorgradient = jet(max(output.exp.sweep));
        subplot(4,4,1:8)
        hold on; grid on
        plot(output.data.time,output.Imon.raw(ctr,:),'color',colorgradient(ctr,:))
    end
    
    subplot(4,4,1:8)
    ylabel('Current [nA]','FontSize',12)
    xlabel('Time [s]','FontSize',12)
    title('VC-Step Protocol','FontSize',22)
    colormap(jet(max(output.exp.sweep)));
    c = colorbar;
    caxis([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(max(output.exp.sweep))])
    c.Label.String = 'Voltage [mV]';
    
    subplot(4,4,9)
    hold on; grid on
    plot(output.data.plots.protocol.stimulusrange,output.data.plots.events,colorplot{1})
    ylabel('Event Count','FontSize',12)
    title('Event count','FontSize',12)
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,10)
    hold on; grid on
    plot(output.data.plots.protocol.stimulusrange,output.data.plots.pks(:,1) + output.data.plots.pksinitiation(:,1),colorplot{1})
    ylabel('Current [nA]','FontSize',12)
    title('1st event amplitude','FontSize',12)
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
  
    subplot(4,4,11)
    hold on; grid on
    plot(output.data.plots.protocol.stimulusrange,output.data.plots.pks(:,2) + output.data.plots.pksinitiation(:,2),colorplot{1})
    ylabel('Current [nA]','FontSize',12)
    title('2nd event amplitude','FontSize',12)
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,12)
    hold on; grid on
    for ctr = 1:size(output.data.plots.pks,1)
        plot(output.data.plots.protocol.stimulusrange(ctr),100 - ((output.data.plots.pks(ctr,2) + output.data.plots.pksinitiation(ctr,2)) / (output.data.plots.pks(ctr,1) + output.data.plots.pksinitiation(ctr,1)) * 100),colorplot{1})
    end
    ylabel('Amp change (%)','FontSize',12)
    title('1-2 amplitude difference %','FontSize',12)
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,13)
    hold on; grid on
    plot(output.data.plots.protocol.stimulusrange,output.data.plots.halfwidth,colorplot{1})
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('half width','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,14)
    hold on; grid on
    for i = 1:size(output.data.plots.pks,1)
        if isnan(output.data.plots.pksindex(i,1)) == 0
            plot(output.data.plots.protocol.stimulusrange(i,1),output.data.time(output.data.plots.pksindex(i,1)-output.data.plots.protocol.stimbeginloc),colorplot{1})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('1st event latency','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,15)
    hold on; grid on
    for i = 1:size(output.data.plots.pks,1)
        if isnan(output.data.plots.pksindex(i,2)) == 0
            plot(output.data.plots.protocol.stimulusrange(i,1),output.data.time(output.data.plots.pksindex(i,2)-output.data.plots.protocol.stimbeginloc),colorplot{1})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('2nd event latency','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
    subplot(4,4,16)
    hold on; grid on
    for i = 1:size(output.data.plots.pks,1)
        if isnan(output.data.plots.pksindex(i,2)) == 0
            plot(output.data.plots.protocol.stimulusrange(i,1),output.data.time(output.data.plots.pksindex(i,2))-output.data.time(output.data.plots.pksindex(i,1)),colorplot{1})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('1-2 latency difference','FontSize',12)
    xlim([output.data.plots.protocol.stimulusrange(1) output.data.plots.protocol.stimulusrange(end)])
    
elseif strcmp(protocol,'VCai') == 1
    
    for ctr = 1:max(output.exp.sweep)
        subplot(5,4,1:4)
        hold on; grid on
        plot(output.data.time(1:2000),output.Vmon.raw(ctr,1:2000),'b')
        plot(output.data.time(2000:3200),output.Vmon.raw(ctr,2000:3200),'r')
    end
    title('VC-Step Protocol: Step (Blue) vs Constant (Red)','FontSize',22)
    ylabel('Voltage [mV]')
    
    for ctr = 1:max(output.exp.sweep)
        colorgradient = jet(max(output.exp.sweep));
        subplot(5,4,5:12)
        hold on; grid on
        plot(output.data.time,output.Imon.raw(ctr,:),'color',colorgradient(ctr,:))
    end
    
    subplot(5,4,5:12)
    ylabel('Current [nA]','FontSize',12)
    xlabel('Time [s]','FontSize',12)
    colormap(jet(max(output.exp.sweep)));
    c = colorbar;
    caxis([output.data.plots.stimulusrange(1) output.data.plots.stimulusrange(max(output.exp.sweep))])
    c.Label.String = 'Voltage [mV]';
        
    subplot(5,4,13)
    hold on; grid on
    plot(output.data.plots.stimulusrange,output.data.plots.step.events,colorplot{1})
    plot(output.data.plots.stimulusrange,output.data.plots.constant.events,colorplot{2})
    ylabel('Event Count','FontSize',12)
    title('Event count','FontSize',12)
    
    subplot(5,4,14)
    hold on; grid on
    plot(output.data.plots.stimulusrange,abs(output.data.plots.step.pks(:,1) + output.data.plots.step.pksinitiation(:,1)),colorplot{1})
    plot(output.data.plots.stimulusrange,abs(output.data.plots.constant.pks(:,1) + output.data.plots.constant.pksinitiation(:,1)),colorplot{2})
    ylabel('Current [nA]','FontSize',12)
    title('1st event amplitude','FontSize',12)
    
    subplot(5,4,15)
    hold on; grid on
    plot(output.data.plots.stimulusrange,output.data.plots.step.pks(:,2) + output.data.plots.step.pksinitiation(:,2),colorplot{1})
    plot(output.data.plots.stimulusrange,output.data.plots.constant.pks(:,2) + output.data.plots.constant.pksinitiation(:,2),colorplot{2})
    ylabel('Current [nA]','FontSize',12)
    title('2nd event amplitude','FontSize',12)
    
    subplot(5,4,16)
    hold on; grid on
    for ctr = 1:size(output.data.plots.step.pks,1)
        plot(output.data.plots.stimulusrange(ctr),100 - ((output.data.plots.step.pks(ctr,2) + output.data.plots.step.pksinitiation(ctr,2)) ./ (output.data.plots.step.pks(ctr,1) + output.data.plots.step.pksinitiation(ctr,1)) * 100),colorplot{1})
    end
    for ctr = 1:size(output.data.plots.constant.pks,1)
        plot(output.data.plots.stimulusrange(ctr),100 - ((output.data.plots.constant.pks(ctr,2) + output.data.plots.constant.pksinitiation(ctr,2)) ./ (output.data.plots.constant.pks(ctr,1) + output.data.plots.constant.pksinitiation(ctr,1)) * 100),colorplot{2})
    end
    ylabel('Amp change (%)','FontSize',12)
    title('1-2 amplitude difference %','FontSize',12)
    
    subplot(5,4,17)
    hold on; grid on    
    plot(output.data.plots.stimulusrange,output.data.plots.step.halfwidth,colorplot{1})
    plot(output.data.plots.stimulusrange,output.data.plots.constant.halfwidth,colorplot{2})
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('half width','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
   
    subplot(5,4,18)
    hold on; grid on
    for i = 1:size(output.data.plots.step.pks,1)
        if isnan(output.data.plots.step.pksindex(i,1)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.step.pksindex(i,1)-output.data.plots.step.stimbeginloc),colorplot{1})
        end
        if isnan(output.data.plots.constant.pksindex(i,1)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.constant.pksindex(i,1)-output.data.plots.constant.stimbeginloc),colorplot{2})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('1st event latency','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
    
    subplot(5,4,19)
    hold on; grid on
    for i = 1:size(output.data.plots.step.pks,1)
        if isnan(output.data.plots.step.pksindex(i,2)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.step.pksindex(i,2)-output.data.plots.step.stimbeginloc),colorplot{1})
        end
        if isnan(output.data.plots.constant.pksindex(i,2)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.constant.pksindex(i,2)-output.data.plots.constant.stimbeginloc),colorplot{2})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('2nd event latency','FontSize',12)
    A = gca;
    A.YAxis.Exponent = 0;
    
    subplot(5,4,20)
    hold on; grid on
    for i = 1:size(output.data.plots.step.pks,1)
        if isnan(output.data.plots.step.pksindex(i,2)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.step.pksindex(i,2))-output.data.time(output.data.plots.step.pksindex(i,1)),colorplot{1})
        end
        if isnan(output.data.plots.constant.pksindex(i,2)) == 0
            plot(output.data.plots.stimulusrange(i,1),output.data.time(output.data.plots.constant.pksindex(i,2))-output.data.time(output.data.plots.constant.pksindex(i,1)),colorplot{2})
        end
    end
    ylabel('Time [s]','FontSize',12)
    xlabel('Membrane potential [mV]','FontSize',12)
    title('1-2 latency difference','FontSize',12)
end

set(gcf,'Position', get(0, 'Screensize'));
if exist([cd filesep 'fig_' datafolder],'dir') == 0
    mkdir(['fig_' datafolder]); 
    savefig(['fig_' datafolder '/' filename])
    saveas(gcf,['fig_' datafolder '/' filename],'png')
else
    savefig(['fig_' datafolder '/' filename])
    saveas(gcf,['fig_' datafolder '/' filename],'png')
end
close all
