function analyze_N2A (filename,savefolder)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_ST_N2A is a function used for a voltage clamp sawtooth protocol
% analysis (N2A version)

% This function retrieves the data exported from HEKA patch clamp setups
% EPC 9 and EPC 10 double. It puts the data into matlab and performs data
% analysis. That data that is used has to be collected during voltage clamp
% experiments.

% Command entry example: analyze_ST_N2A('170926_AB_ST100_C')
close all; 

load ([filename '.mat']);

files = whos ('Trace*'); % lists the current variables, in long format form
output.file.name = filename;
disp (['analyze_ST_N2A: ' output.file.name]) % displays what you are doing, the VC_Step_Analysis

%% Organization of HEKA patch clamp data:

% Trace_X_Y_Z_W
% X: experiment number
% Y: protocol number
% Z: sweep
% W: channels.  

% Example: 
% Experiment AB8 with the protocols CC_C,CC_C,CC_C,ST100_C,ST50_C,ST10_C,VC_C in that order
% analyze_ST_N2A('170926_AB_ST100_C')
% Resulting in Trace_X_Y_Z_W
% X = 8     (AB8)
% Y = 4     (ST100_C is in position 4 in the chain of protocols shown in line 30
% Z = 1-2/3  (ST100 protocol consist of 2/3 repetitions, Z will change from 1 to 3)
% W = 1-4   (Setup 1 has 4 channels (Imon2-Vmon2-Imon1-Vmon1), setup 2 has
%           2 channels (Vmon1-Imon1). As this experiment is performed on setup 1, we
%           will have 4 different channels. W will change from 1 to 4.
% Resulting in Trace_8_4_1_1 for the first sweep in channel 1

%% Channel configuration

% The channel configuration depends on which setup u use. As setup 1 has
% two amplifiers it will give four output channels (Vmon1 Vmon2 Imon1 Imon2.
% Setup 2 and 3 have both one amplifier so the channel configuration is
% always: column 1 = Imon, column 2 = Vmon. When using amplifier 2 on setup
% one u use the samp configuration. Only when using amplifier 1 on setup 1
% you need to use: column 3 = Imon, column 4 = Vmon

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

for ctr = 1:3
    output.data.sweep{ctr}.time {1} = output.data.values{1}(:,1); % as the time is always the same, it does not matter which trace we pick
end

%% Stimulation configuration

% Stimbegin depends on the protocol that is being used. For now, only the ST100 is used. Next versions 
% will also have the stimbegin of the 200ms, 50ms, 10ms sawtooths.
stimbegin = 0.05:0.2:0.85;      % ST protocols consist of 5 sawtooths starting at 0.05, 0.25, 0.45, 0.65, and 0.85 seconds.
mindistance = 0.19;             % mindistance is used to help the first findpeaks function to find the main spikes.

% Sweepstartindex is calculated with the help of stimbegin which changes
% for each protocol. Stimbegin is shown in line 7.
for i = 1:3
    for ctr = 1:5
        output.data.sweep{i}.plots.sweepstartindex(ctr) = find(abs(output.data.sweep{i}.time{1}-stimbegin(ctr)) == min(abs(output.data.sweep{i}.time{1}-stimbegin(ctr))));
        output.data.sweep{i}.plots.sweepstart(ctr) = output.data.sweep{i}.time{1}(output.data.sweep{i}.plots.sweepstartindex(ctr));
    end
end

%% Split and reorder the data.sweep{ctr}. 
 
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

%% Prelocations

for ctr = 1:max(output.exp.sweep)
    % Output.data.sweep{ctr}.plots is used for the main spikes 
    % of each sawtooth (the first sodium channel spike). These are used to set
    % some parameters and calculate
    output.data.sweep{ctr}.plots.pksindex = zeros(1,length(stimbegin));
    output.data.sweep{ctr}.plots.pksinitiation = zeros(1,length(stimbegin));
    output.data.sweep{ctr}.plots.locAfterStim = zeros(1,length(stimbegin));
    output.data.sweep{ctr}.plots.halfwidth = zeros(1,length(stimbegin));

    % output.data.sweep{ctr}.findpeaks is used for all the other peaks. This is done
    % seperatly so findpeaks can search for smaller peaks between big events.
    % As the protocol is repeated 2 or 3 times, the struct is divided into 3
    % different sweeps.
    output.data.sweep{ctr}.findpeaks.sweep = zeros(5,10);                   
    output.data.sweep{ctr}.findpeaks.sweeploc = zeros(5,10);    
    output.data.sweep{ctr}.findpeaks.sweepprom = nan(5,10);       
    output.data.sweep{ctr}.findpeaks.sweepwidth = zeros(5,10);    
    output.data.sweep{ctr}.findpeaks.sweepperdiff = nan(5,10); 
end

%% Collecting/analysis of the dat
for ctr = 1:max(output.exp.sweep)   % loops over the analysis for the amount of sweeps (2 or 3).
 
    % The function findpeaks from the signalling toolbox is used to find
    % the main spikes of each sawotooths. This is done by implementing a
    % minimal distance between peaks. The mindistance is a bit smaller than
    % the step size of the stimulus.
    
    for tooth = 1:5
        [output.data.sweep{ctr}.plots.pks(1,tooth),output.data.sweep{ctr}.plots.loc(1,tooth)] = findpeaks(-output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.sweepstartindex(tooth) + 500):(output.data.sweep{ctr}.plots.sweepstartindex(tooth) + 4000)),output.data.sweep{ctr}.time{1}((output.data.sweep{ctr}.plots.sweepstartindex(tooth) + 500):(output.data.sweep{ctr}.plots.sweepstartindex(tooth) + 4000)),'MinPeakProminence',0.025,'MinPeakDistance',0.05);
    end
    
    [output.data.sweep{ctr}.plots.pkssaw(1,:),output.data.sweep{ctr}.plots.locsaw(1,:)] = findpeaks(output.Imon.raw(ctr,:),output.data.sweep{ctr}.time{1},'MinPeakDistance',mindistance);

    % Colleting of important parameters like the initiationtime of the
    % peak, index of pks and initiation, location of the spike relative to
    % the stimulus.
    for i = 1:length(output.data.sweep{ctr}.plots.pks)
        output.data.sweep{ctr}.plots.pksindex(1,i)= find(-output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex(i)+500:end) == output.data.sweep{ctr}.plots.pks(1,i),1) + output.data.sweep{ctr}.plots.sweepstartindex(i) + 500 - 1;
        output.data.sweep{ctr}.plots.pkssawindex(1,i)= find(output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.pksindex(1,i):end)) == output.data.sweep{ctr}.plots.pkssaw(1,i),1) + output.data.sweep{ctr}.plots.pksindex(1,i)- 1;
        output.data.sweep{ctr}.plots.pksinitiation(1,i)= max(output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.sweepstartindex(i)+500):output.data.sweep{ctr}.plots.pksindex(1,i)));
        output.data.sweep{ctr}.plots.pksinitiationindex(1,i)= find(output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex(i):output.data.sweep{ctr}.plots.pksindex(1,i)) == output.data.sweep{ctr}.plots.pksinitiation(1,i),1) + output.data.sweep{ctr}.plots.sweepstartindex(i);
        output.data.sweep{ctr}.plots.locAfterStim(1,i)= output.data.sweep{ctr}.plots.loc(1,i)- (mean(diff(stimbegin)) * (i-1) + stimbegin(1));
       
        % Calculation of the halfwidth of the main spike.
        hh = abs((-output.Imon.raw(ctr,:)) - (output.data.sweep{ctr}.plots.pks(1,i)- ((output.data.sweep{ctr}.plots.pks(1,i)+ output.data.sweep{ctr}.plots.pksinitiation(1,i)) / 2)));
        output.data.sweep{ctr}.plots.halfwidth1(1,i)= find(hh((output.data.sweep{ctr}.plots.pksinitiationindex(1,i)):output.data.sweep{ctr}.plots.pksindex(1,i)) == min(hh((output.data.sweep{ctr}.plots.pksinitiationindex(1,i)):output.data.sweep{ctr}.plots.pksindex(1,i)),[],2),1) + output.data.sweep{ctr}.plots.pksinitiationindex(1,i);
        output.data.sweep{ctr}.plots.halfwidth2(1,i)= (output.data.sweep{ctr}.plots.halfwidth1(1,i)- 1) + (find(output.Imon.raw(ctr,output.data.sweep{ctr}.plots.halfwidth1(1,i) + 50:end) > output.Imon.raw(ctr,output.data.sweep{ctr}.plots.halfwidth1(1,i)),1));
    end
    
    % Now all the data around every spike in the measurment is collected.
    % An if/else construction is used to place data from sweep 1 in
    % output.data.sweep{ctr}.findpeaks.sweep etc..
    for peak = 1:length(output.data.sweep{ctr}.plots.pks)
        [intE,intEl,intEw,intEp] = findpeaks(-output.Imon.raw(ctr,output.data.sweep{ctr}.plots.pksinitiationindex(1,peak):output.data.sweep{ctr}.plots.pkssawindex(1,peak)),'MinPeakProminence',0.05,'MinPeakDistance',100);
        output.data.sweep{ctr}.findpeaks.sweep(peak,1:size(intE,2)) = intE;
        output.data.sweep{ctr}.findpeaks.sweeploc(peak,1:size(intEl,2)) = intEl  + output.data.sweep{ctr}.plots.pksinitiationindex(1,peak) - 1;
        output.data.sweep{ctr}.findpeaks.sweepwidth(peak,1:size(intE,2)) = intEw;
        output.data.sweep{ctr}.findpeaks.sweepprom(peak,1:size(intE,2)) = intEp;
        output.data.sweep{ctr}.findpeaks.sweepevents(1,peak) = nnz(output.data.sweep{ctr}.findpeaks.sweeploc(peak,:));
        if output.data.sweep{ctr}.findpeaks.sweeploc(peak,2) == 0
            output.data.sweep{ctr}.findpeaks.sweeptimediff(1,peak) = NaN;
        else
            output.data.sweep{ctr}.findpeaks.sweeptimediff(1,peak) = output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.findpeaks.sweeploc(peak,2)) - output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.findpeaks.sweeploc(peak,1));
        end
    end
    
    % Special loop for the percentage difference between the first and
    % second peak.
    if isnan(output.data.sweep{ctr}.findpeaks.sweepprom(1,1)) == 0
        for tooth = 1:5
            for peak = 1:4
                output.data.sweep{ctr}.findpeaks.sweepperdiff(tooth,peak) = 100 - (output.data.sweep{ctr}.findpeaks.sweepprom(tooth,peak+1) / output.data.sweep{ctr}.findpeaks.sweepprom(tooth,peak) * 100);
            end
        end
    else    
        output.data.sweep{ctr}.findpeaks.sweepperdiff = zeros(5,4);
    end
    % If there are no spikes (good example are the TTX measurments)
    % there is nothing to calculate. The struct will theb be filled
    % with NaN's. Note: some parts of the struct are not mentioned here
    % as those were already prelocated with NaN's or zeros.
    if isempty(output.data.sweep{ctr}.plots.pks)
        output.data.sweep{ctr}.plots.pks(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.loc(1,1:5) = NaN;
        output.data.sweep{ctr}.findpeaks.sweeptimediff(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.pksinitiationindex(1,1:5) = NaN;
        output.data.sweep{ctr}.findpeaks.sweepevents(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.currentVStime(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.halfwidth(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.initiation(1,1:5) = NaN;
    else
        output.data.sweep{ctr}.plots.initiation(1,:) = output.Vmon.raw(ctr,output.data.sweep{ctr}.plots.pksinitiationindex(1,:));
        output.data.sweep{ctr}.plots.halfwidth(1,:) = output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.halfwidth2(1,:)) - output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.halfwidth1(1,:));
        output.data.sweep{ctr}.plots.currentVStime(1,:) = (output.data.sweep{ctr}.plots.pks(1,:) + output.data.sweep{ctr}.plots.pksinitiation(1,:))./((output.data.sweep{ctr}.plots.loc(1,:)-(output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.pksinitiationindex(1,:)))')*1000);
    end
end

if exist(savefolder,'dir') == 0
    mkdir(savefolder); save([savefolder '/analyzed_ST_' filename], 'output')        
else
    save([savefolder '/analyzed_ST_' filename], 'output')
end
