function analyze_VCstep(filename,savefolder)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_VCstep is part of the mother function 'analyze_patch' and is a 
% function used for a voltage clamp sawtooth protocol analysis

% This function retrieves the data exported from HEKA patch clamp setups
% EPC 9 and EPC 10 double. It puts the data into matlab and performs data
% analysis. That data that is used has to be collected during voltage clamp
% experiments.

% Command entry example: analyze_VCstep('KK_MagnetoSlices_12012018_E1_5_VClampStep_WithoutMagnet')

close all; 

load ([filename '.mat']);

files = whos ('Trace*');                        % lists the current variables, in long format form
output.file.name = filename;
disp (['analyze_VCstep: ' output.file.name])    % displays what you are doing, the VC_Step_Analysis

%% Organization of HEKA patch clamp data:

% Trace_X_Y_Z_W
% X: experiment number
% Y: protocol number
% Z: sweep
% W: channels.  

% Example: 
% Experiment AB8 with the protocols CC_C,CC_C,CC_C,ST100_C,ST50_C,ST10_C,VC_C in that order
% analyze_VCstep('170926_AB_VCstep')
% Resulting in Trace_X_Y_Z_W
% X = 8     (AB8)
% Y = 7     (VC_C is in position 7 in the chain)
% Z = 1-2/3 (Setup 1 has 4 channels (Imon2-Vmon2-Imon1-Vmon1), setup 2 has
% W = 1-4   2 channels (Vmon1-Imon1). As this experiment is performed on setup 1, we
%           will have 4 different channels. W will change from 1 to 4.
% Resulting in Trace_8_7_1_1 for the first sweep in channel 1

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
    output.data.prevalues {Tc} = eval(files(Tc).name);
    output.data.sweepduration (Tc) =  output.data.prevalues {Tc}(end,1)-output.data.prevalues{Tc}(1,1);             % in seconds
    output.data.samplingrate (Tc) =  round(size(output.data.prevalues {Tc},1)/output.data.sweepduration (Tc));   % in Hz
end

output.data.time(:,1) = output.data.prevalues{1}(:,1); % as the time is always the same, it does not matter which trace we pick

%% Trace ordering
% This reordering is needed because matlab is special........ not the good kind.

[~,order] = sort(output.exp.sweep);
output.data.values = output.data.prevalues(order);

%% Channel configuration

% Depending on the setup, the traces will consist of multiple channels.
% Here we let the code know where to look for the current and voltge data.sweep{ctr}.
if round(max(output.data.values{2}(:,2)),2) == 0    
    Ichannel = 3;
    Vchannel = 4;
else
    Ichannel = 1;
    Vchannel = 2;
end

%% Split and reorder the data. 
 
% First we will put all the Vmon data (Column 2 from each sweep) together.
% The same is done for the Imon data (Column 1 from each sweep. The
% position in which it is placed is changed by the counters(Vctr / Ictr) as
% after each placement the counter ads 1, so the position moves from n to
% n+1. 
 
Vctr = 1; Ictr = 1; % Vctr and Ictr function as counter

for Tc = 1:size(files,1) 
    
    if output.exp.channel(Tc) == Vchannel % only picks from column 2, as we assigned 2 = Vchannel
        output.Vmon.raw(Vctr,:) = 1000.*output.data.values{Tc}(:,2); % Data is V, we want it in mV so we multiply with 1000.
    
        Vctr = Vctr + 1;
           
    elseif output.exp.channel(Tc) == Ichannel % only picks from column 1, as we assigned 1 = Vchannel
        output.Imon.raw(Ictr,:) = 10e8.*output.data.values{Tc}(:,2); % Data is A, we want it in nA so we multiply with 10e8.
    
        Ictr = Ictr + 1;
    end
end

%% Stimulation configuration

% Stimbegin depends on the protocol that is being used. This code only
% supports VC_C

output.data.plots.protocol.stimbeginloc = find(diff(output.Vmon.raw(1,:)) == max(diff(output.Vmon.raw(1,:)))) - 1;
output.data.plots.protocol.stimbegin = output.data.time(output.data.plots.protocol.stimbeginloc);
output.data.plots.protocol.stimendloc = find(diff(output.Vmon.raw(1,:)) == min(diff(output.Vmon.raw(1,:)))) - 1;
output.data.plots.protocol.stimend = output.data.time(output.data.plots.protocol.stimendloc);
output.data.plots.protocol.stimulusrange = round(output.Vmon.raw(:,output.data.plots.protocol.stimbeginloc+25)); % output.data.plots.stimulusrange is the membrane potentials in the protocol. Stimbegin + 25 is used as the amplifier needs a few ms to get to the target level.
output.data.plots.protocol.sweepstartindex = find(abs(output.data.time-output.data.plots.protocol.stimbegin) == min(abs(output.data.time-output.data.plots.protocol.stimbegin)));
output.data.plots.protocol.sweepstart = output.data.time(output.data.plots.protocol.sweepstartindex);

clearvars Tc Ictr Vctr % clean up some of the non-important variables.

%% Collecting/analysis of the data

% Prelocation of the spike peaks(pks) and location(loc). Maximum of 25
% peaks for each sweep is placed to limit the size. When there are more
% then 25 spikes an error will occur. Just make the number higher then.
output.data.plots.pks = nan(size(output.Vmon.raw,1),25);
output.data.plots.loc = nan(size(output.Vmon.raw,1),25);
output.data.plots.pksindex = nan(size(output.Vmon.raw,1),25);
output.data.plots.pksinitiation = nan(size(output.Vmon.raw,1),25);
output.data.plots.pksinitiationindex = nan(size(output.Vmon.raw,1),25);

for ctr = 1:(size(output.Vmon.raw,1))

% Findpeaks is finding the main (biggest starting peaks) peaks in the
% sawtooth. The minimum peak height is set on 15% of the biggest
% amplitude change in the VCstep. If not all the main peaks are
% found, lower the MinPeakHeight.     
 
    [pks,loc] = findpeaks(-output.Imon.raw(ctr,1:output.data.plots.protocol.stimendloc),output.data.time(1:output.data.plots.protocol.stimendloc),'MinpeakProminence',0.5);
    
    if length(pks) > 1
        if pks(1) < pks(2)
            pks = []; loc = [];
            [pks,loc] =  findpeaks(-output.Imon.raw(ctr,output.data.plots.protocol.sweepstartindex + 20:output.data.plots.protocol.stimendloc),output.data.time(output.data.plots.protocol.sweepstartindex + 20:output.data.plots.protocol.stimendloc),'MinpeakProminence',0.15);
        end
    end
    
    output.data.plots.pks(ctr,1:length(pks)) = pks;
    output.data.plots.loc(ctr,1:length(loc)) = loc';
    output.data.plots.events(ctr,1) = sum(~isnan(output.data.plots.pks(ctr,:)));

    for i = 1:size(output.data.plots.pks,2)
        if isnan(output.data.plots.pks(ctr,i)) == 0
            h = find(-output.Imon.raw(ctr,output.data.plots.protocol.sweepstartindex:end) == output.data.plots.pks(ctr,i)) + output.data.plots.protocol.sweepstartindex - 1;
            while i > 1 && h(1) <= output.data.plots.pksindex(ctr,i-1)
                h(1) = [];
            end
            output.data.plots.pksindex(ctr,i) = h(1);
            
            if i == 1
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.protocol.sweepstartindex:output.data.plots.protocol.sweepstartindex + 20) == max(output.Imon.raw(ctr,output.data.plots.protocol.sweepstartindex:output.data.plots.protocol.sweepstartindex + 20)),1) + output.data.plots.protocol.sweepstartindex - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.pksindex(ctr,i)),output.data.plots.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
            
            if i == 2
                stimpeakindex = find(output.Imon.raw(ctr,output.data.plots.pksindex(ctr,1):output.data.plots.pksindex(ctr,2)) == max(output.Imon.raw(ctr,output.data.plots.pksindex(ctr,1):output.data.plots.pksindex(ctr,2))),1) + output.data.plots.pksindex(ctr,1) - 1;
                initiationtrace = output.Imon.raw(ctr,stimpeakindex:output.data.plots.pksindex(ctr,i));
                initiationline = linspace(output.Imon.raw(ctr,stimpeakindex),output.Imon.raw(ctr,output.data.plots.pksindex(ctr,i)),output.data.plots.pksindex(ctr,i)-stimpeakindex+1);
                leftovers = initiationtrace-initiationline;
                output.data.plots.pksinitiationindex(ctr,i) = find(leftovers == max(leftovers),1) + stimpeakindex - 1;
                output.data.plots.pksinitiation(ctr,i) = output.Imon.raw(ctr,output.data.plots.pksinitiationindex(ctr,i));
                clearvars stimpeakindex initiationtrace initiationline leftovers
            end
        end
    end
end

for i = 1:(size(output.Vmon.raw,1))     
    if isnan(output.data.plots.pks(i,1))
        output.data.plots.halfwidth(1,i) = NaN;
    else
        hh = abs((-output.Imon.raw(i,:)) - (output.data.plots.pks(i,1)- ((output.data.plots.pks(i,1)+ output.data.plots.pksinitiation(i,1)) / 2)));
        output.data.plots.halfwidth1(1,i)= find(hh(output.data.plots.pksinitiationindex(i,1):output.data.plots.pksindex(i,1)) == min(hh(output.data.plots.pksinitiationindex(i,1):output.data.plots.pksindex(i,1)),[],2)) + output.data.plots.pksinitiationindex(i,1);
        output.data.plots.halfwidth2(1,i)= (output.data.plots.halfwidth1(1,i)- 1) + (find(output.Imon.raw(i,output.data.plots.halfwidth1(1,i):end) > output.Imon.raw(i,output.data.plots.halfwidth1(1,i)),1));
        output.data.plots.halfwidth(1,i) = output.data.time(output.data.plots.halfwidth2(1,i)) - output.data.time(output.data.plots.halfwidth1(1,i));
      end
end

if exist([cd filesep savefolder],'dir') == 0
    mkdir(savefolder); save([savefolder '/analyzed_VCstep_' filename], 'output')        
else
    save([savefolder '/analyzed_VCstep_' filename], 'output')
end