function analyze_VCai(filename,savefolder)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_AI is part of the mother function 'analyze_patch' and is a
% function used for voltage clamp activat/inactivate protocol analysis.

% This function retrieves the data exported from HEKA patch clamp setups
% EPC 9 and EPC 10 double. It puts the data into matlab and performs data
% analysis. That data that is used has to be collected during voltage clamp
% experiments.

% Command entry example: analyze_AI('MP_acute_slices_ethanol_control_030418_p53_ate_ActInact')

close all

load ([filename '.mat']);

files = whos ('Trace*');                    % lists the current variables, in long format form
output.file.name = filename;
disp (['analyze_AI: ' output.file.name])    % displays what you are doing, the VC_Step_Analysis

%% Organization of HEKA patch clamp data:

% Trace_X_Y_Z_W
% X: experiment number
% Y: protocol number
% Z: sweep
% W: channels.  

% Example: 
% Experiment AB8 with the protocols CC_C,CC_C,CC_C,ST100_C,ST50_C,ST10_C,AI in that order
% analyze_AI('170926_AB_AI')
% Resulting in Trace_X_Y_Z_W
% X = 8     (AB8)
% Y = 7     (AI is in position position 7 in the chain)
% Z = 1-2/3 (Setup 1 has 4 channels (Imon2-Vmon2-Imon1-Vmon1), setup 2 has
% W = 1-4   2 channels (Vmon1-Imon1). As this experiment is performed on setup 1, we
%           will have 4 different channels. W will change from 1 to 4.
% Resulting in Trace_8_7_1_1 for the first sweep in channel 1

%% Channel configuration

% Depending on the setup, the traces will consist of multiple channels.
% Here we let the code know where to look for the current and voltge data.sweep{ctr}.
Ichannel = 1;
Vchannel = 2;

%% Organizing the data

% Tc = TraceCount, Ts = TraceSplit,  ctr = counter

for Tc = 1:size(files,1) % loops until every trace is used, amount of traces is equal to the size(files,1). As it lists all the traces.
    
    Ts = strsplit(files(Tc).name,'_'); % Splits the traces at the '_'. For our example it will split in 'Trace', '4', '7', '1', '1'.
    
% Organizing and Assigning the split traces. str2double is used to convert
% the elements in the cell array of character to vectors C to double    
    
    output.exp.experiment (Tc) = str2double(Ts{2});
    output.exp.dataset (Tc) = str2double(Ts{3});
    output.exp.sweep (Tc) = str2double(Ts{4});
    output.exp.channel (Tc) = str2double(Ts{5});

    output.data.tracename {Tc} = files(Tc).name;
    output.data.values {Tc} = eval(files(Tc).name);
    output.data.sweepduration (Tc) =  output.data.values {Tc}(end,1)-output.data.values{Tc}(1,1);             % in seconds
    output.data.samplingrate (Tc) =  round(size(output.data.values {Tc},1)/output.data.sweepduration (Tc));   % in Hz
end

output.data.time(:,1) = output.data.values{1}(:,1); % as the time is always the same, it does not matter which trace we pick

%% Stimulation configuration

% Stimbegin depends on the protocol that is being used. This code only
% supports AI

stimbegin = 0.0500;                              % Stimulation begins at time point 0.050 seconds. 
output.data.plots.step.stimbeginloc = 1000;      % This is index 1000 in output.data.time{1}.
output.data.plots.constant.stimbeginloc = 2000;                                      
stimbeginloc = 1000;
stimendloc = 3000;                               % stimend = 0.150 sec; 
stimdivider = 2000;                              % As the protocol consist of two different parts, the increasing step pulse and the constant hold. The function will look at it seperatly. The divider shows the point the step ends and the constant begins.

%% Split and reorder the data. 
 
% First we will put all the Vmon data (Column 2 from each sweep) together.
% The same is done for the Imon data (Column 1 from each sweep. The
% position in which it is placed is changed by the counters(Vctr / Ictr) as
% after each placement the counter ads 1, so the position moves from n to
% n+1. 

Vctr = 1; Ictr = 1; % Vctr and Ictr function as counter

for Tc = 1:size(files,1) 
    
    if output.exp.channel(Tc) == Vchannel % only picks from column 2, as we assigned 2 = Vchannel
        output.Vmon.pre(Vctr,:) = 1000.*output.data.values{Tc}(:,2); % Data is V, we want it in mV so we multiply with 1000.
    
        Vctr = Vctr + 1;
           
    elseif output.exp.channel(Tc) == Ichannel % only picks from column 1, as we assigned 1 = Vchannel
        output.Imon.pre(Ictr,:) = 10e8.*output.data.values{Tc}(:,2); % Data is A, we want it in nA so we multiply with 10e8.
    
        Ictr = Ictr + 1;
    end
end

% Sorry for this, it is dirty, but it works.
if max(output.exp.sweep) == 12
    output.Vmon.raw(1,:) = output.Vmon.pre(4,:);         output.Imon.raw(1,:) = output.Imon.pre(4,:);
    output.Vmon.raw(2,:) = output.Vmon.pre(5,:);         output.Imon.raw(2,:) = output.Imon.pre(5,:);
    output.Vmon.raw(3,:) = output.Vmon.pre(6,:);         output.Imon.raw(3,:) = output.Imon.pre(6,:);
    output.Vmon.raw(4,:) = output.Vmon.pre(7,:);         output.Imon.raw(4,:) = output.Imon.pre(7,:);
    output.Vmon.raw(5,:) = output.Vmon.pre(8,:);         output.Imon.raw(5,:) = output.Imon.pre(8,:);
    output.Vmon.raw(6,:) = output.Vmon.pre(9,:);         output.Imon.raw(6,:) = output.Imon.pre(9,:);
    output.Vmon.raw(7,:) = output.Vmon.pre(10,:);        output.Imon.raw(7,:) = output.Imon.pre(10,:);
    output.Vmon.raw(8,:) = output.Vmon.pre(11,:);        output.Imon.raw(8,:) = output.Imon.pre(11,:);
    output.Vmon.raw(9,:) = output.Vmon.pre(12,:);        output.Imon.raw(9,:) = output.Imon.pre(12,:);
    output.Vmon.raw(10,:) = output.Vmon.pre(1,:);        output.Imon.raw(10,:) = output.Imon.pre(1,:);
    output.Vmon.raw(11,:) = output.Vmon.pre(2,:);        output.Imon.raw(11,:) = output.Imon.pre(2,:);
    output.Vmon.raw(12,:) = output.Vmon.pre(3,:);        output.Imon.raw(12,:) = output.Imon.pre(3,:);
elseif max(output.exp.sweep) == 13
    output.Vmon.raw(1,:) = output.Vmon.pre(5,:);         output.Imon.raw(1,:) = output.Imon.pre(5,:);
    output.Vmon.raw(2,:) = output.Vmon.pre(6,:);         output.Imon.raw(2,:) = output.Imon.pre(6,:);
    output.Vmon.raw(3,:) = output.Vmon.pre(7,:);         output.Imon.raw(3,:) = output.Imon.pre(7,:);
    output.Vmon.raw(4,:) = output.Vmon.pre(8,:);         output.Imon.raw(4,:) = output.Imon.pre(8,:);
    output.Vmon.raw(5,:) = output.Vmon.pre(9,:);         output.Imon.raw(5,:) = output.Imon.pre(9,:);
    output.Vmon.raw(6,:) = output.Vmon.pre(10,:);        output.Imon.raw(6,:) = output.Imon.pre(10,:);
    output.Vmon.raw(7,:) = output.Vmon.pre(11,:);        output.Imon.raw(7,:) = output.Imon.pre(11,:);
    output.Vmon.raw(8,:) = output.Vmon.pre(12,:);        output.Imon.raw(8,:) = output.Imon.pre(12,:);
    output.Vmon.raw(9,:) = output.Vmon.pre(13,:);        output.Imon.raw(9,:) = output.Imon.pre(13,:);
    output.Vmon.raw(10,:) = output.Vmon.pre(1,:);        output.Imon.raw(10,:) = output.Imon.pre(1,:);
    output.Vmon.raw(11,:) = output.Vmon.pre(2,:);        output.Imon.raw(11,:) = output.Imon.pre(2,:);
    output.Vmon.raw(12,:) = output.Vmon.pre(3,:);        output.Imon.raw(12,:) = output.Imon.pre(3,:);
    output.Vmon.raw(13,:) = output.Vmon.pre(4,:);        output.Imon.raw(13,:) = output.Imon.pre(4,:);
elseif max(output.exp.sweep) == 14
    output.Vmon.raw(1,:) = output.Vmon.pre(6,:);         output.Imon.raw(1,:) = output.Imon.pre(6,:);
    output.Vmon.raw(2,:) = output.Vmon.pre(7,:);         output.Imon.raw(2,:) = output.Imon.pre(7,:);
    output.Vmon.raw(3,:) = output.Vmon.pre(8,:);         output.Imon.raw(3,:) = output.Imon.pre(8,:);
    output.Vmon.raw(4,:) = output.Vmon.pre(9,:);         output.Imon.raw(4,:) = output.Imon.pre(9,:);
    output.Vmon.raw(5,:) = output.Vmon.pre(10,:);        output.Imon.raw(5,:) = output.Imon.pre(10,:);
    output.Vmon.raw(6,:) = output.Vmon.pre(11,:);        output.Imon.raw(6,:) = output.Imon.pre(11,:);
    output.Vmon.raw(7,:) = output.Vmon.pre(12,:);        output.Imon.raw(7,:) = output.Imon.pre(12,:);
    output.Vmon.raw(8,:) = output.Vmon.pre(13,:);        output.Imon.raw(8,:) = output.Imon.pre(13,:);
    output.Vmon.raw(9,:) = output.Vmon.pre(14,:);        output.Imon.raw(9,:) = output.Imon.pre(14,:);
    output.Vmon.raw(10,:) = output.Vmon.pre(1,:);        output.Imon.raw(10,:) = output.Imon.pre(1,:);
    output.Vmon.raw(11,:) = output.Vmon.pre(2,:);        output.Imon.raw(11,:) = output.Imon.pre(2,:);
    output.Vmon.raw(12,:) = output.Vmon.pre(3,:);        output.Imon.raw(12,:) = output.Imon.pre(3,:);
    output.Vmon.raw(13,:) = output.Vmon.pre(4,:);        output.Imon.raw(13,:) = output.Imon.pre(4,:);
    output.Vmon.raw(14,:) = output.Vmon.pre(5,:);        output.Imon.raw(14,:) = output.Imon.pre(5,:);    
end

output.data.plots.stimulusrange = round(output.Vmon.raw(:,output.data.plots.step.stimbeginloc+25)); % output.data.plots.stimulusrange is the membrane potentials in the protocol. Stimbegin + 25 is used as the amplifier needs a few ms to get to the target level.

clearvars Tc Ictr Vctr % clean up some of the non-important variables.

%% Collecting/analysis of the data

% Sweepstartindex is calculated with the help of stimbegin which changes
% for each protocol. 
output.data.plots.sweepstartindex = find(abs(output.data.time-stimbegin) == min(abs(output.data.time-stimbegin)));
output.data.plots.sweepstart = output.data.time(output.data.plots.sweepstartindex);

% Prelocation of the spike peaks(pks) and location(loc). Maximum of 25
% peaks for each sweep is placed to limit the size. When there are more
% then 25 spikes an error will occur. Just make the number highet then.
output.data.plots.step.pks = nan(size(output.Vmon.raw,1),25);
output.data.plots.step.loc = nan(size(output.Vmon.raw,1),25);
output.data.plots.step.pksindex = nan(size(output.Vmon.raw,1),25);
output.data.plots.step.pksinitiation = nan(size(output.Vmon.raw,1),25);
output.data.plots.step.pksinitiationindex = nan(size(output.Vmon.raw,1),25);

output.data.plots.constant.pks = nan(size(output.Vmon.raw,1),25);
output.data.plots.constant.loc = nan(size(output.Vmon.raw,1),25);
output.data.plots.constant.pksindex = nan(size(output.Vmon.raw,1),25);
output.data.plots.constant.pksinitiation = nan(size(output.Vmon.raw,1),25);
output.data.plots.constant.pksinitiationindex = nan(size(output.Vmon.raw,1),25);

for ctr = 1:(size(output.Vmon.raw,1))

    % Findpeaks is finding the main (biggest starting peaks) peaks in the
    % sawtooth. The minimum peak height is set on 15% of the biggest
    % amplitude change in the VCstep. If not all the main peaks are
    % found, lower the MinPeakHeight. 
    
    [pksS,locS] = findpeaks(-output.Imon.raw(ctr,1:stimdivider),output.data.time(1:stimdivider),'MinpeakProminence',0.5);
    
    [pksC,locC] = findpeaks(-output.Imon.raw(ctr,stimdivider:stimendloc),output.data.time(stimdivider:stimendloc),'MinpeakProminence',0.5);

    if length(pksS) > 1
        if pksS(1) < pksS(2)
            pksS = []; locS = [];
            [pksS,locS] =  findpeaks(-output.Imon.raw(ctr,output.data.plots.sweepstartindex + 20:stimdivider),output.data.time(output.data.plots.sweepstartindex + 20:stimdivider),'MinpeakProminence',0.15);
        end
    end
    
    if length(pksC) > 1
        if pksC(1) < pksC(2)
            pksC = []; locC = [];
            [pksC,locC] =  findpeaks(-output.Imon.raw(ctr,stimdivider:stimendloc),output.data.time(stimdivider:stimendloc),'MinpeakProminence',0.05);
        end
    end
    
    output.data.plots.step.pks(ctr,1:length(pksS)) = pksS;
    output.data.plots.step.loc(ctr,1:length(locS)) = locS';
    output.data.plots.step.events(ctr,1) = sum(~isnan(output.data.plots.step.pks(ctr,:)));
    output.data.plots.constant.pks(ctr,1:length(pksC)) = pksC;
    output.data.plots.constant.loc(ctr,1:length(locC)) = locC';
    output.data.plots.constant.events(ctr,1) = sum(~isnan(output.data.plots.constant.pks(ctr,:)));
    
    for i = 1:size(output.data.plots.step.pks,2)
        if isnan(output.data.plots.step.pks(ctr,i)) == 0
            h = find(-output.Imon.raw(ctr,output.data.plots.sweepstartindex:end) == output.data.plots.step.pks(ctr,i)) + output.data.plots.sweepstartindex - 1;
            while i > 1 && h(1) <= output.data.plots.step.pksindex(ctr,i-1)
                h(1) = [];
            end
            output.data.plots.step.pksindex(ctr,i) = h(1);
            
            if i == 1
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.sweepstartindex:output.data.plots.sweepstartindex + 20) == max(output.Imon.raw(ctr,output.data.plots.sweepstartindex:output.data.plots.sweepstartindex + 20)),1) + output.data.plots.sweepstartindex - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.step.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.step.pksindex(ctr,i)),output.data.plots.step.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.step.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.step.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.step.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
            
            if i == 2
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.step.pksindex(ctr,1):output.data.plots.step.pksindex(ctr,2)) == max(output.Imon.raw(ctr,output.data.plots.step.pksindex(ctr,1):output.data.plots.step.pksindex(ctr,2))),1) + output.data.plots.step.pksindex(ctr,1) - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.step.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.step.pksindex(ctr,i)),output.data.plots.step.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.step.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.step.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.step.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
        end
    end
    
    for i = 1:size(output.data.plots.constant.pks,2)
        if isnan(output.data.plots.constant.pks(ctr,i)) == 0
            h = find(-output.Imon.raw(ctr,output.data.plots.sweepstartindex:end) == output.data.plots.constant.pks(ctr,i)) + output.data.plots.sweepstartindex - 1;
            while i > 1 && h(1) <= output.data.plots.constant.pksindex(ctr,i-1)
                h(1) = [];
            end
            output.data.plots.constant.pksindex(ctr,i) = h(1);
            
            if i == 1
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.sweepstartindex:output.data.plots.sweepstartindex + 20) == max(output.Imon.raw(ctr,output.data.plots.sweepstartindex:output.data.plots.sweepstartindex + 20)),1) + output.data.plots.sweepstartindex - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.constant.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.constant.pksindex(ctr,i)),output.data.plots.constant.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.constant.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.constant.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.constant.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
            
            if i == 2
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.constant.pksindex(ctr,1):output.data.plots.constant.pksindex(ctr,2)) == max(output.Imon.raw(ctr,output.data.plots.constant.pksindex(ctr,1):output.data.plots.constant.pksindex(ctr,2))),1) + output.data.plots.constant.pksindex(ctr,1) - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.constant.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.constant.pksindex(ctr,i)),output.data.plots.constant.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.constant.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.constant.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.constant.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
        end
    end
end

for i = 1:(size(output.Vmon.raw,1))     
    if isnan(output.data.plots.step.pks(i,1))
        output.data.plots.step.halfwidth(1,i) = NaN;
    else
        hh = abs((-output.Imon.raw(i,:)) - (output.data.plots.step.pks(i,1)- ((output.data.plots.step.pks(i,1)+ output.data.plots.step.pksinitiation(i,1)) / 2)));
        output.data.plots.step.halfwidth1(1,i)= find(hh(output.data.plots.step.pksinitiationindex(i,1):output.data.plots.step.pksindex(i,1)) == min(hh(output.data.plots.step.pksinitiationindex(i,1):output.data.plots.step.pksindex(i,1)),[],2)) + output.data.plots.step.pksinitiationindex(i,1);
        output.data.plots.step.halfwidth2(1,i)= (output.data.plots.step.halfwidth1(1,i)- 1) + (find(output.Imon.raw(i,output.data.plots.step.halfwidth1(1,i):end) > output.Imon.raw(i,output.data.plots.step.halfwidth1(1,i)),1));
        output.data.plots.step.halfwidth(1,i) = output.data.time(output.data.plots.step.halfwidth2(1,i)) - output.data.time(output.data.plots.step.halfwidth1(1,i));
    end
end

for i = 1:(size(output.Vmon.raw,1))     
    if isnan(output.data.plots.constant.pks(i,1))
        output.data.plots.constant.halfwidth(1,i) = NaN;
    else
        hh = abs((-output.Imon.raw(i,:)) - (output.data.plots.constant.pks(i,1)- ((output.data.plots.constant.pks(i,1)+ output.data.plots.constant.pksinitiation(i,1)) / 2)));
        output.data.plots.constant.halfwidth1(1,i)= find(hh(output.data.plots.constant.pksinitiationindex(i,1):output.data.plots.constant.pksindex(i,1)) == min(hh(output.data.plots.constant.pksinitiationindex(i,1):output.data.plots.constant.pksindex(i,1)),[],2)) + output.data.plots.constant.pksinitiationindex(i,1);
        output.data.plots.constant.halfwidth2(1,i)= (output.data.plots.constant.halfwidth1(1,i)- 1) + (find(output.Imon.raw(i,output.data.plots.constant.halfwidth1(1,i):end) > output.Imon.raw(i,output.data.plots.constant.halfwidth1(1,i)),1));
        output.data.plots.constant.halfwidth(1,i) = output.data.time(output.data.plots.constant.halfwidth2(1,i)) - output.data.time(output.data.plots.constant.halfwidth1(1,i));
    end
end

if exist([cd filesep savefolder],'dir') == 0
    mkdir(savefolder); save([savefolder '/analyzed_VCai_' filename], 'output')        
else
    save([savefolder '/analyzed_VCai_' filename], 'output')
end

