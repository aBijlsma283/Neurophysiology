function analyze_Celine (group1,group2,group3,group4,group5,group6)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

%% Setup & Protocol additions
Ichannel = 1;  
Vchannel = 2;
stimbegin = 0.05;

%% Grouping of data
filesgroup1 = dir (['Celine*' (group1) '.mat']);
filesgroup2 = dir (['Celine*' (group2) '.mat']);
filesgroup3 = dir (['Celine*' (group3) '.mat']);
filesgroup4 = dir (['Celine*' (group4) '.mat']);
filesgroup5 = dir (['Celine*' (group5) '.mat']);
filesgroup6 = dir (['Celine*' (group6) '.mat']);
filesgroup7 = dir (['EJKK*' (group1) '.mat']);
filesgroup8 = dir (['EJKK*' (group2) '.mat']);
filesgroup9 = dir (['EJKK*' (group3) '.mat']);
filesgroup10= dir (['EJKK*' (group4) '.mat']);
filesgroup11= dir (['EJKK*' (group5) '.mat']);
filesgroup12= dir (['EJKK*' (group6) '.mat']);

lijn1 = '_10_';
lijn2 = '_23_';

filesgroup13= dir (['Celine*' lijn1 '*ST50(80mV).mat']);
filesgroup14= dir (['Celine*' lijn1 '*ST50(100mV).mat']);
filesgroup15= dir (['Celine*' lijn1 '*ST100(80mV).mat']);
filesgroup16= dir (['Celine*' lijn1 '*ST100(100mV).mat']);
filesgroup17= dir (['Celine*' lijn1 '*ST200(80mV).mat']);
filesgroup18= dir (['Celine*' lijn1 '*ST200(100mV).mat']);
filesgroup19= dir (['Celine*' lijn2 '*ST50(80mV).mat']);
filesgroup20= dir (['Celine*' lijn2 '*ST50(100mV).mat']);
filesgroup21= dir (['Celine*' lijn2 '*ST100(80mV).mat']);
filesgroup22= dir (['Celine*' lijn2 '*ST100(100mV).mat']);
filesgroup23= dir (['Celine*' lijn2 '*ST200(80mV).mat']);
filesgroup24= dir (['Celine*' lijn2 '*ST200(100mV).mat']);

%% Failure lists
Failure_Celine_Errorlist = {};
Failure_WT_Errorlist = {};

%% Data seperation & analysis
for GroupCount = 1:24
    
    clearvars output
    
    if GroupCount == 1 || GroupCount == 2 || GroupCount == 7 || GroupCount == 8 || GroupCount == 13 || GroupCount == 14 || GroupCount == 19 || GroupCount == 20
       stimhelp = 250;
       stimend = 2000;
    elseif GroupCount == 3 || GroupCount == 4 || GroupCount == 9 || GroupCount == 10 || GroupCount == 15 || GroupCount == 16 || GroupCount == 21 || GroupCount == 22
       stimhelp = 500;
       stimend = 2500;
    else
       stimhelp = 1000;
       stimend = 3500;
    end    
    
    gdata{GroupCount}.pks(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    gdata{GroupCount}.loc(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    gdata{GroupCount}.pksindex(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    gdata{GroupCount}.pksinitiation(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    gdata{GroupCount}.pksinitiationindex(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    gdata{GroupCount}.pksdepolarizationindex(size(eval(['filesgroup' num2str(GroupCount)]),1),3) = 0;
    
    for FileCount = 1:size(eval(['filesgroup' num2str(GroupCount)]),1)
        name = strsplit(eval(['filesgroup' num2str(GroupCount) '(' num2str(FileCount) ')' '.name']));
        load (name{1});
        
        files = whos ('Trace*');
        
        for Count3 = 1:size(files,1)
            
            Ts = strsplit(files(Count3).name,'_');
            
            output.exp.experiment (Count3) = str2double(Ts{2});
            output.exp.dataset (Count3) = str2double(Ts{3});
            output.exp.sweep (Count3) = str2double(Ts{4});
            output.exp.channel (Count3) = str2double(Ts{5});
            
            output.data.tracename {Count3} = files(Count3).name;
            output.data.values {Count3} = eval(files(Count3).name);
            output.data.sweepduration (Count3) =  output.data.values {Count3}(end,1)-output.data.values{Count3}(1,1);
            output.data.samplingrate (Count3) =  round(size(output.data.values {Count3},1)/output.data.sweepduration (Count3));
        end
        
        for i = 1:3
            output.data.sweep{i}.time {1} = output.data.values{1}(:,1);
            output.data.sweep{i}.plots.sweepstartindex(1) = find(abs(output.data.sweep{i}.time{1}-stimbegin) == min(abs(output.data.sweep{i}.time{1}-stimbegin)));
        end
        
        try
            CountVctr = 1; CountIctr = 1;
            
            for Count4 = 1:size(files,1)
                
                if output.exp.channel(Count4) == Vchannel
                    output.Vmon.raw(CountVctr,:) = 1000.*output.data.values{Count4}(:,2);
                    
                    CountVctr = CountVctr + 1;
                    
                elseif output.exp.channel(Count4) == Ichannel
                    output.Imon.raw(CountIctr,:) = 10e8.*output.data.values{Count4}(:,2);
                    
                    CountIctr = CountIctr + 1;
                end
            end
            
            for ctr = 1:max(output.exp.sweep)
                [gdata{GroupCount}.pks(FileCount,ctr),gdata{GroupCount}.loc(FileCount,ctr)] = findpeaks(-output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.sweepstartindex(1) + stimhelp):stimend),output.data.sweep{ctr}.time{1}((output.data.sweep{ctr}.plots.sweepstartindex(1) + stimhelp):stimend),'MinPeakProminence',0.025);
                gdata{GroupCount}.pksindex(FileCount,ctr)= find(-output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex + stimhelp:stimend) == gdata{GroupCount}.pks(FileCount,ctr),1) + output.data.sweep{ctr}.plots.sweepstartindex + stimhelp - 1;
                gdata{GroupCount}.pksinitiation(FileCount,ctr)= max(output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.sweepstartindex + stimhelp):gdata{GroupCount}.pksindex(FileCount,ctr)));
                gdata{GroupCount}.pksinitiationindex(FileCount,ctr)= find(output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex:gdata{GroupCount}.pksindex(FileCount,ctr)) == gdata{GroupCount}.pksinitiation(FileCount,ctr),1) + output.data.sweep{ctr}.plots.sweepstartindex - 1;
                gdata{GroupCount}.prominence = gdata{GroupCount}.pks + gdata{GroupCount}.pksinitiation;
                if gdata{GroupCount}.pksinitiationindex(FileCount,ctr) ~= 0
                    gdata{GroupCount}.pksdepolarizationindex(FileCount,ctr) = find(-output.Imon.raw(ctr,gdata{GroupCount}.pksindex(FileCount,ctr):end) <= -gdata{GroupCount}.pksinitiation(FileCount,ctr),1) + gdata{GroupCount}.pksindex(FileCount,ctr) - 1;
                end
            end
        catch
            if GroupCount <= 6
                Failure_Celine_Errorlist{end+1} = name{1};
            else
                Failure_WT_Errorlist{end+1} = name{1};
            end
        end
                       
        clear -regexp ^Trace_
    end
    
    gdata{GroupCount}.initiationtimeindex = gdata{GroupCount}.pksindex - gdata{GroupCount}.pksinitiationindex;
    gdata{GroupCount}.initiationtimeindex(find(gdata{GroupCount}.initiationtimeindex == 0)) = 1;
    gdata{GroupCount}.initiationtime = output.data.sweep{1}.time{1}(gdata{GroupCount}.initiationtimeindex);
    output.Vmon.initiation = output.Vmon.raw(1,:)';
    gdata{GroupCount}.pksinitiationindex(find(gdata{GroupCount}.pksinitiationindex == 0)) = 1;
    gdata{GroupCount}.initiationvoltage = output.Vmon.initiation(gdata{GroupCount}.pksinitiationindex);
    gdata{GroupCount}.initiationvoltage(find(gdata{GroupCount}.pksinitiationindex == 1)) = NaN;
    gdata{GroupCount}.pksdepolarizationindex(find(gdata{GroupCount}.pksdepolarizationindex == 0)) = 2;
    gdata{GroupCount}.depolarizationindex = gdata{GroupCount}.pksdepolarizationindex - gdata{GroupCount}.pksinitiationindex;
    gdata{GroupCount}.depolarizationindex(find(gdata{GroupCount}.depolarizationindex == 0)) = 1;
    gdata{GroupCount}.depolarizationtime = output.data.sweep{1}.time{1}(gdata{GroupCount}.depolarizationindex);
end
keyboard
%% Statistic analysis
[gdata{25}.promH(1),gdata{25}.promP(1)] = ttest2(gdata{1}.prominence(:),gdata{7}.prominence(:));
[gdata{25}.promH(2),gdata{25}.promP(2)] = ttest2(gdata{2}.prominence(:),gdata{8}.prominence(:));
[gdata{25}.promH(3),gdata{25}.promP(3)] = ttest2(gdata{3}.prominence(:),gdata{9}.prominence(:));
[gdata{25}.promH(4),gdata{25}.promP(4)] = ttest2(gdata{4}.prominence(:),gdata{10}.prominence(:));
[gdata{25}.promH(5),gdata{25}.promP(5)] = ttest2(gdata{5}.prominence(:),gdata{11}.prominence(:));
[gdata{25}.promH(6),gdata{25}.promP(6)] = ttest2(gdata{6}.prominence(:),gdata{12}.prominence(:));

[gdata{25}.locH(1),gdata{25}.locP(1)] = ttest2(gdata{1}.loc(:),gdata{7}.loc(:));
[gdata{25}.locH(2),gdata{25}.locP(2)] = ttest2(gdata{2}.loc(:),gdata{8}.loc(:));
[gdata{25}.locH(3),gdata{25}.locP(3)] = ttest2(gdata{3}.loc(:),gdata{9}.loc(:));
[gdata{25}.locH(4),gdata{25}.locP(4)] = ttest2(gdata{4}.loc(:),gdata{10}.loc(:));
[gdata{25}.locH(5),gdata{25}.locP(5)] = ttest2(gdata{5}.loc(:),gdata{11}.loc(:));
[gdata{25}.locH(6),gdata{25}.locP(6)] = ttest2(gdata{6}.loc(:),gdata{12}.loc(:));

[gdata{25}.iniTimeH(1),gdata{25}.iniTimeP(1)] = ttest2(gdata{1}.initiationtime(:),gdata{7}.initiationtime(:));
[gdata{25}.iniTimeH(2),gdata{25}.iniTimeP(2)] = ttest2(gdata{2}.initiationtime(:),gdata{8}.initiationtime(:));
[gdata{25}.iniTimeH(3),gdata{25}.iniTimeP(3)] = ttest2(gdata{3}.initiationtime(:),gdata{9}.initiationtime(:));
[gdata{25}.iniTimeH(4),gdata{25}.iniTimeP(4)] = ttest2(gdata{4}.initiationtime(:),gdata{10}.initiationtime(:));
[gdata{25}.iniTimeH(5),gdata{25}.iniTimeP(5)] = ttest2(gdata{5}.initiationtime(:),gdata{11}.initiationtime(:));
[gdata{25}.iniTimeH(6),gdata{25}.iniTimeP(6)] = ttest2(gdata{6}.initiationtime(:),gdata{12}.initiationtime(:));

[gdata{25}.iniVoltH(1),gdata{25}.iniVoltP(1)] = ttest2(gdata{1}.initiationvoltage(:),gdata{7}.initiationvoltage(:));
[gdata{25}.iniVoltH(2),gdata{25}.iniVoltP(2)] = ttest2(gdata{2}.initiationvoltage(:),gdata{8}.initiationvoltage(:));
[gdata{25}.iniVoltH(3),gdata{25}.iniVoltP(3)] = ttest2(gdata{3}.initiationvoltage(:),gdata{9}.initiationvoltage(:));
[gdata{25}.iniVoltH(4),gdata{25}.iniVoltP(4)] = ttest2(gdata{4}.initiationvoltage(:),gdata{10}.initiationvoltage(:));
[gdata{25}.iniVoltH(5),gdata{25}.iniVoltP(5)] = ttest2(gdata{5}.initiationvoltage(:),gdata{11}.initiationvoltage(:));
[gdata{25}.iniVoltH(6),gdata{25}.iniVoltP(6)] = ttest2(gdata{6}.initiationvoltage(:),gdata{12}.initiationvoltage(:));

[gdata{25}.DepoH(1),gdata{25}.DepoP(1)] = ttest2(gdata{1}.depolarizationtime(:),gdata{7}.depolarizationtime(:));
[gdata{25}.DepoH(2),gdata{25}.DepoP(2)] = ttest2(gdata{2}.depolarizationtime(:),gdata{8}.depolarizationtime(:));
[gdata{25}.DepoH(3),gdata{25}.DepoP(3)] = ttest2(gdata{3}.depolarizationtime(:),gdata{9}.depolarizationtime(:));
[gdata{25}.DepoH(4),gdata{25}.DepoP(4)] = ttest2(gdata{4}.depolarizationtime(:),gdata{10}.depolarizationtime(:));
[gdata{25}.DepoH(5),gdata{25}.DepoP(5)] = ttest2(gdata{5}.depolarizationtime(:),gdata{11}.depolarizationtime(:));
[gdata{25}.DepoH(6),gdata{25}.DepoP(6)] = ttest2(gdata{6}.depolarizationtime(:),gdata{12}.depolarizationtime(:));

[gdata{25}.OVpromH(1),gdata{25}.OVpromP(1)] = ttest2(gdata{13}.prominence(:),gdata{19}.prominence(:));
[gdata{25}.OVpromH(2),gdata{25}.OVpromP(2)] = ttest2(gdata{14}.prominence(:),gdata{20}.prominence(:));
[gdata{25}.OVpromH(3),gdata{25}.OVpromP(3)] = ttest2(gdata{15}.prominence(:),gdata{21}.prominence(:));
[gdata{25}.OVpromH(4),gdata{25}.OVpromP(4)] = ttest2(gdata{16}.prominence(:),gdata{22}.prominence(:));
[gdata{25}.OVpromH(5),gdata{25}.OVpromP(5)] = ttest2(gdata{17}.prominence(:),gdata{23}.prominence(:));
[gdata{25}.OVpromH(6),gdata{25}.OVpromP(6)] = ttest2(gdata{18}.prominence(:),gdata{24}.prominence(:));
save('analyzed_N2A_Thesis', 'gdata', 'Failure_Celine_Errorlist', 'Failure_WT_Errorlist')

%% Plotting (General figures)
figure(1)
colorarray = {'ro','bo','go','mo','ko','co'};
subplot(2,2,1)
for plotter = 1:6
    plot(gdata{plotter}.loc,gdata{plotter}.prominence,colorarray{plotter});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{group1,group2,group3,group4,group5,group6})

subplot(2,2,2)
for plotter = 1:6
    plot(gdata{plotter}.initiationtime,gdata{plotter}.prominence,colorarray{plotter});
    hold on
end
xlabel('Duration of depolarization [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{group1,group2,group3,group4,group5,group6})

subplot(2,2,3)
for plotter = 1:6
    plot(gdata{plotter}.loc,gdata{plotter}.initiationvoltage,colorarray{plotter});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{group1,group2,group3,group4,group5,group6})

subplot(2,2,4)
for plotter = 1:6
    plot(gdata{plotter}.initiationtime,gdata{plotter}.initiationvoltage,colorarray{plotter});
    hold on
end
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{group1,group2,group3,group4,group5,group6})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, General Overexpression analysis.png')
savefig('N2A, General Overexpression analysis')
close all

figure(2)
colorarray = {'ro','bo','go','mo','ko','co'};
subplot(2,2,1)
for plotter = 7:12
    plot(gdata{plotter}.loc,gdata{plotter}.prominence,colorarray{plotter-6});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{[group1 '-WT'],[group2 '-WT'],[group3 '-WT'],[group4 '-WT'],[group5 '-WT'],[group6 '-WT']})

subplot(2,2,2)
for plotter = 7:12
    plot(gdata{plotter}.initiationtime,gdata{plotter}.prominence,colorarray{plotter-6});
    hold on
end
xlabel('Duration of depolarization [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{[group1 '-WT'],[group2 '-WT'],[group3 '-WT'],[group4 '-WT'],[group5 '-WT'],[group6 '-WT']})

subplot(2,2,3)
for plotter = 7:12
    plot(gdata{plotter}.loc,gdata{plotter}.initiationvoltage,colorarray{plotter-6});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{[group1 '-WT'],[group2 '-WT'],[group3 '-WT'],[group4 '-WT'],[group5 '-WT'],[group6 '-WT']})

subplot(2,2,4)
for plotter = 7:12
    plot(gdata{plotter}.initiationtime,gdata{plotter}.initiationvoltage,colorarray{plotter-6});
    hold on
end
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); h3 = plot(NaN,NaN,'go'); h4 = plot(NaN,NaN,'mo'); h5 = plot(NaN,NaN,'ko'); h6 = plot(NaN,NaN,'co'); 
legend([h1 h2 h3 h4 h5 h6],{[group1 '-WT'],[group2 '-WT'],[group3 '-WT'],[group4 '-WT'],[group5 '-WT'],[group6 '-WT']})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, General WT analysis.png')
savefig('N2A, General WT analysis')
close all

%% Plotting (parameter vs parameter)
figure(3)
subplot(1,3,1)
plot(gdata{1}.loc,gdata{1}.prominence,colorarray{1});
hold on
plot(gdata{2}.loc,gdata{2}.prominence,colorarray{2});
plot(gdata{7}.loc,gdata{7}.prominence,colorarray{3});
plot(gdata{8}.loc,gdata{8}.prominence,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('50ms Amplitude vs Timing')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group1,group2,[group1 '-WT'],[group2 '-WT']})

subplot(1,3,2)
plot(gdata{3}.loc,gdata{3}.prominence,colorarray{1});
hold on
plot(gdata{4}.loc,gdata{4}.prominence,colorarray{2});
plot(gdata{9}.loc,gdata{9}.prominence,colorarray{3});
plot(gdata{10}.loc,gdata{10}.prominence,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('100ms Amplitude vs Timing')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group3,group4,[group3 '-WT'],[group4 '-WT']})

subplot(1,3,3)
plot(gdata{5}.loc,gdata{5}.prominence,colorarray{1});
hold on
plot(gdata{6}.loc,gdata{6}.prominence,colorarray{2});
plot(gdata{11}.loc,gdata{11}.prominence,colorarray{3});
plot(gdata{12}.loc,gdata{12}.prominence,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('200ms Amplitude vs Timing')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group5,group6,[group5 '-WT'],[group6 '-WT']})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Amplitude vs Timing.png')
savefig('N2A, Amplitude vs Timing')
close all

figure(4)
subplot(1,3,1)
plot(gdata{1}.initiationtime,gdata{1}.initiationvoltage,colorarray{1});
hold on
plot(gdata{2}.initiationtime,gdata{2}.initiationvoltage,colorarray{2});
plot(gdata{7}.initiationtime,gdata{7}.initiationvoltage,colorarray{3});
plot(gdata{8}.initiationtime,gdata{8}.initiationvoltage,colorarray{4});
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('50ms Initiation voltage vs Slope duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group1,group2,[group1 '-WT'],[group2 '-WT']})

subplot(1,3,2)
plot(gdata{3}.initiationtime,gdata{3}.initiationvoltage,colorarray{1});
hold on
plot(gdata{4}.initiationtime,gdata{4}.initiationvoltage,colorarray{2});
plot(gdata{9}.initiationtime,gdata{9}.initiationvoltage,colorarray{3});
plot(gdata{10}.initiationtime,gdata{10}.initiationvoltage,colorarray{4});
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('100ms Initiation voltage vs Slope duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group3,group4,[group3 '-WT'],[group4 '-WT']})

subplot(1,3,3)
plot(gdata{5}.initiationtime,gdata{5}.initiationvoltage,colorarray{1});
hold on
plot(gdata{6}.initiationtime,gdata{6}.initiationvoltage,colorarray{2});
plot(gdata{11}.initiationtime,gdata{11}.initiationvoltage,colorarray{3});
plot(gdata{12}.initiationtime,gdata{12}.initiationvoltage,colorarray{4});
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('200ms Initiation voltage vs Slope duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group5,group6,[group5 '-WT'],[group6 '-WT']})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Initiation voltage vs Duration till peak.png')
savefig('N2A, Initiation voltage vs Duration till peak')
close all

figure(5)
subplot(1,3,1)
plot(gdata{1}.loc,gdata{1}.depolarizationtime,colorarray{1});
hold on
plot(gdata{2}.loc,gdata{2}.depolarizationtime,colorarray{2});
plot(gdata{7}.loc,gdata{7}.depolarizationtime,colorarray{3});
plot(gdata{8}.loc,gdata{8}.depolarizationtime,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Depolarization duration [mS]','FontSize',16)
title('50ms Timing vs Depolarization duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group1,group2,[group1 '-WT'],[group2 '-WT']})

subplot(1,3,2)
plot(gdata{3}.loc,gdata{3}.depolarizationtime,colorarray{1});
hold on
plot(gdata{4}.loc,gdata{4}.depolarizationtime,colorarray{2});
plot(gdata{9}.loc,gdata{9}.depolarizationtime,colorarray{3});
plot(gdata{10}.loc,gdata{10}.depolarizationtime,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Depolarization duration [mS]','FontSize',16)
title('100ms Timing vs Depolarization duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group3,group4,[group3 '-WT'],[group4 '-WT']})

subplot(1,3,3)
plot(gdata{5}.loc,gdata{5}.depolarizationtime,colorarray{1});
hold on
plot(gdata{6}.loc,gdata{6}.depolarizationtime,colorarray{2});
plot(gdata{11}.loc,gdata{11}.depolarizationtime,colorarray{3});
plot(gdata{12}.loc,gdata{12}.depolarizationtime,colorarray{4});
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Depolarization duration [mS]','FontSize',16)
title('200ms Timing vs Depolarization duration')

a1 = plot(NaN,NaN,'ro'); a2 = plot(NaN,NaN,'bo'); a3 = plot(NaN,NaN,'go'); a4 = plot(NaN,NaN,'mo'); 
legend([a1 a2 a3 a4],{group5,group6,[group5 '-WT'],[group6 '-WT']})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Timing vs Depolarization duration.png')
savefig('N2A, Timing vs Depolarization duration')
close all

figure(6)
colorarray2 = {'ro','ro','ro','ro','ro','ro','bo','bo','bo','bo','bo','bo'};
subplot(2,2,1)
for plotter = 13:24
    plot(gdata{plotter}.loc,gdata{plotter}.prominence,colorarray2{plotter-12});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); 
legend([h1 h2],{'10','23'})

subplot(2,2,2)
for plotter = 13:24
    plot(gdata{plotter}.initiationtime,gdata{plotter}.prominence,colorarray2{plotter-12});
    hold on
end
xlabel('Duration of depolarization [mS]','FontSize',16)
ylabel('Current amplitude [nA]','FontSize',16)
title('Amplitude vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo');
legend([h1 h2],{'10','23'})

subplot(2,2,3)
for plotter = 13:24
    plot(gdata{plotter}.loc,gdata{plotter}.initiationvoltage,colorarray2{plotter-12});
    hold on
end
xlabel('Time after Stimuli [mS]','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Timing')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); 
legend([h1 h2],{'10','23'})

subplot(2,2,4)
for plotter = 13:24
    plot(gdata{plotter}.initiationtime,gdata{plotter}.initiationvoltage,colorarray2{plotter-12});
    hold on
end
xlabel('Duration of depolarization','FontSize',16)
ylabel('Initiation voltage [mV]','FontSize',16)
title('Initiation voltage vs Depolarization duration')
h1 = plot(NaN,NaN,'ro'); h2 = plot(NaN,NaN,'bo'); 
legend([h1 h2],{'10','23'})

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Overexpression 10 vs 23.png')
savefig('N2A, Overexpression 10 vs 23')
close all

%% Plotting (Boxplots)
for ctr = 1:24
    gdata{ctr}.prominence(find(gdata{ctr}.prominence == 0)) = NaN;
    gdata{ctr}.depolarizationtime(find(gdata{ctr}.depolarizationtime == 0)) = NaN;
    gdata{ctr}.initiationvoltage(find(gdata{ctr}.initiationvoltage == 0)) = NaN;
end

marker_size = 6;
figure(7)
b1 = subplot(2,3,1); hold on;
for ctr = 1:size(filesgroup7,1)
    plot(1,gdata{7}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup1,1)
    plot(2,gdata{1}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

b2 = subplot(2,3,2); hold on;
for ctr = 1:size(filesgroup9,1)
    plot(1,gdata{9}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup3,1)
    plot(2,gdata{3}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

b3 = subplot(2,3,3); hold on;
for ctr = 1:size(filesgroup11,1)
    plot(1,gdata{11}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup5,1)
    plot(2,gdata{5}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

b4 = subplot(2,3,4); hold on;
for ctr = 1:size(filesgroup8,1)
    plot(1,gdata{8}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup2,1)
    plot(2,gdata{2}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

b5 = subplot(2,3,5); hold on;
for ctr = 1:size(filesgroup10,1)
    plot(1,gdata{10}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup4,1)
    plot(2,gdata{4}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

b6 = subplot(2,3,6); hold on;
for ctr = 1:size(filesgroup12,1)
    plot(1,gdata{12}.prominence(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup6,1)
    plot(2,gdata{6}.prominence(ctr,1),'ro','MarkerSize', marker_size)
end

boxplot_position1 = ones(50,1)*0.5;
boxplot_position2 = ones(50,1)*2.5;
boxplot_width = [0.3, 0.3];

axes(b1); hold on
boxplot(gdata{7}.prominence(:,1) , boxplot_position1(1:size(gdata{7}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{1}.prominence(:,1) , boxplot_position2(1:size(gdata{1}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group1 ' p = ' num2str(gdata{25}.promP(1,1))])

axes(b2); hold on
boxplot(gdata{9}.prominence(:,1) , boxplot_position1(1:size(gdata{9}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{3}.prominence(:,1) , boxplot_position2(1:size(gdata{3}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group3 ' p = ' num2str(gdata{25}.promP(1,3))])

axes(b3); hold on
boxplot(gdata{11}.prominence(:,1) , boxplot_position1(1:size(gdata{11}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{5}.prominence(:,1) , boxplot_position2(1:size(gdata{5}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group5 ' p = ' num2str(gdata{25}.promP(1,5))])

axes(b4); hold on
boxplot(gdata{8}.prominence(:,1) , boxplot_position1(1:size(gdata{8}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{2}.prominence(:,1) , boxplot_position2(1:size(gdata{2}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group2 ' p = ' num2str(gdata{25}.promP(1,2))])

axes(b5); hold on
boxplot(gdata{10}.prominence(:,1) , boxplot_position1(1:size(gdata{10}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{4}.prominence(:,1) , boxplot_position2(1:size(gdata{4}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group4 ' p = ' num2str(gdata{25}.promP(1,4))])

axes(b6); hold on
boxplot(gdata{12}.prominence(:,1) , boxplot_position1(1:size(gdata{12}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{6}.prominence(:,1) , boxplot_position2(1:size(gdata{6}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group6 ' p = ' num2str(gdata{25}.promP(1,6))])

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Boxplots Prominence.png')
savefig('N2A, Boxplots Prominence')
close all

figure(8)
b1 = subplot(2,3,1); hold on;
for ctr = 1:size(filesgroup7,1)
    plot(1,gdata{7}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup1,1)
    plot(2,gdata{1}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

b2 = subplot(2,3,2); hold on;
for ctr = 1:size(filesgroup9,1)
    plot(1,gdata{9}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup3,1)
    plot(2,gdata{3}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

b3 = subplot(2,3,3); hold on;
for ctr = 1:size(filesgroup11,1)
    plot(1,gdata{11}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup5,1)
    plot(2,gdata{5}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

b4 = subplot(2,3,4); hold on;
for ctr = 1:size(filesgroup8,1)
    plot(1,gdata{8}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup2,1)
    plot(2,gdata{2}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

b5 = subplot(2,3,5); hold on;
for ctr = 1:size(filesgroup10,1)
    plot(1,gdata{10}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup4,1)
    plot(2,gdata{4}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

b6 = subplot(2,3,6); hold on;
for ctr = 1:size(filesgroup12,1)
    plot(1,gdata{12}.depolarizationtime(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup6,1)
    plot(2,gdata{6}.depolarizationtime(ctr,1),'ro','MarkerSize', marker_size)
end

axes(b1); hold on
boxplot(gdata{7}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{7}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{1}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{1}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group1 ' p = ' num2str(gdata{25}.DepoP(1,1))])

axes(b2); hold on
boxplot(gdata{9}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{9}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{3}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{3}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group3 ' p = ' num2str(gdata{25}.DepoP(1,3))])

axes(b3); hold on
boxplot(gdata{11}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{11}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{5}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{5}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group5 ' p = ' num2str(gdata{25}.DepoP(1,5))])

axes(b4); hold on
boxplot(gdata{8}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{8}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{2}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{2}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group2 ' p = ' num2str(gdata{25}.DepoP(1,2))])

axes(b5); hold on
boxplot(gdata{10}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{10}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{4}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{4}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group4 ' p = ' num2str(gdata{25}.DepoP(1,4))])

axes(b6); hold on
boxplot(gdata{12}.depolarizationtime(:,1) , boxplot_position1(1:size(gdata{12}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{6}.depolarizationtime(:,1) , boxplot_position2(1:size(gdata{6}.depolarizationtime,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group6 ' p = ' num2str(gdata{25}.DepoP(1,6))])

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Boxplots depolarizationtime.png')
savefig('N2A, Boxplots depolarizationtime')
close all

figure(9)
b1 = subplot(2,3,1); hold on;
for ctr = 1:size(filesgroup7,1)
    plot(1,gdata{7}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup1,1)
    plot(2,gdata{1}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

b2 = subplot(2,3,2); hold on;
for ctr = 1:size(filesgroup9,1)
    plot(1,gdata{9}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup3,1)
    plot(2,gdata{3}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

b3 = subplot(2,3,3); hold on;
for ctr = 1:size(filesgroup11,1)
    plot(1,gdata{11}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup5,1)
    plot(2,gdata{5}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

b4 = subplot(2,3,4); hold on;
for ctr = 1:size(filesgroup8,1)
    plot(1,gdata{8}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup2,1)
    plot(2,gdata{2}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

b5 = subplot(2,3,5); hold on;
for ctr = 1:size(filesgroup10,1)
    plot(1,gdata{10}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup4,1)
    plot(2,gdata{4}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

b6 = subplot(2,3,6); hold on;
for ctr = 1:size(filesgroup12,1)
    plot(1,gdata{12}.initiationvoltage(ctr,1),'bo','MarkerSize', marker_size)
end
for ctr = 1:size(filesgroup6,1)
    plot(2,gdata{6}.initiationvoltage(ctr,1),'ro','MarkerSize', marker_size)
end

axes(b1); hold on
boxplot(gdata{7}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{7}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{1}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{1}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group1 ' p = ' num2str(gdata{25}.iniVoltP(1,1))])

axes(b2); hold on
boxplot(gdata{9}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{9}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{3}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{3}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group3 ' p = ' num2str(gdata{25}.iniVoltP(1,3))])

axes(b3); hold on
boxplot(gdata{11}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{11}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{5}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{5}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group5 ' p = ' num2str(gdata{25}.iniVoltP(1,5))])

axes(b4); hold on
boxplot(gdata{8}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{8}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{2}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{2}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group2 ' p = ' num2str(gdata{25}.iniVoltP(1,2))])

axes(b5); hold on
boxplot(gdata{10}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{10}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{4}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{4}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group4 ' p = ' num2str(gdata{25}.iniVoltP(1,4))])

axes(b6); hold on
boxplot(gdata{12}.initiationvoltage(:,1) , boxplot_position1(1:size(gdata{12}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 0.5 , 'Widths', boxplot_width);
boxplot(gdata{6}.initiationvoltage(:,1) , boxplot_position2(1:size(gdata{6}.initiationvoltage,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2.5 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2],'xticklabel',{'WT','Overexpression'})
ylabel('Current [nA]')
title([group6 ' p = ' num2str(gdata{25}.iniVoltP(1,6))])

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Boxplots initiationvoltage.png')
savefig('N2A, Boxplots initiationvoltage')
close all
keyboard
figure(9)
boxplot_position3 = ones(50,1)*1;
boxplot_position4 = ones(50,1)*2;
boxplot_position5 = ones(50,1)*3;

subplot(2,3,1); hold on
boxplot(gdata{7}.prominence(:,1) , boxplot_position3(1:size(gdata{7}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{13}.prominence(:,1) , boxplot_position4(1:size(gdata{13}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{19}.prominence(:,1) , boxplot_position5(1:size(gdata{19}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group1)

subplot(2,3,2); hold on
boxplot(gdata{9}.prominence(:,1) , boxplot_position3(1:size(gdata{9}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{15}.prominence(:,1) , boxplot_position4(1:size(gdata{15}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{21}.prominence(:,1) , boxplot_position5(1:size(gdata{21}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group3)

subplot(2,3,3); hold on
boxplot(gdata{11}.prominence(:,1) , boxplot_position3(1:size(gdata{11}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{17}.prominence(:,1) , boxplot_position4(1:size(gdata{17}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{23}.prominence(:,1) , boxplot_position5(1:size(gdata{23}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group5)

subplot(2,3,4); hold on
boxplot(gdata{8}.prominence(:,1) , boxplot_position3(1:size(gdata{8}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{14}.prominence(:,1) , boxplot_position4(1:size(gdata{14}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{20}.prominence(:,1) , boxplot_position5(1:size(gdata{20}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group2)

subplot(2,3,5); hold on
boxplot(gdata{10}.prominence(:,1) , boxplot_position3(1:size(gdata{10}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{16}.prominence(:,1) , boxplot_position4(1:size(gdata{16}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{22}.prominence(:,1) , boxplot_position5(1:size(gdata{22}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group4)

subplot(2,3,6); hold on
boxplot(gdata{12}.prominence(:,1) , boxplot_position3(1:size(gdata{12}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{18}.prominence(:,1) , boxplot_position4(1:size(gdata{18}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{24}.prominence(:,1) , boxplot_position5(1:size(gdata{24}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group6)

set(gcf,'Position', get(0, 'Screensize'));
saveas(gcf,'N2A, Boxplots overexpression.png')
savefig('N2A, Boxplots overexpression')
close all

keyboard

figure(10); hold on;
boxplot(gdata{9}.prominence(:,1) , boxplot_position3(1:size(gdata{9}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 1 , 'Widths', boxplot_width);
boxplot(gdata{15}.prominence(:,1) , boxplot_position4(1:size(gdata{15}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 2 , 'Widths', boxplot_width);
boxplot(gdata{21}.prominence(:,1) , boxplot_position5(1:size(gdata{21}.prominence,1),1) , 'MedianStyle', 'line', 'OutlierSize', 4, 'Positions', 3 , 'Widths', boxplot_width);
xlim('auto'); ylim('auto')
set(gca,'xtick',[1 2 3],'xticklabel',{'WT','OE-10','OE-23'})
ylabel('Current [nA]')
title(group3)
