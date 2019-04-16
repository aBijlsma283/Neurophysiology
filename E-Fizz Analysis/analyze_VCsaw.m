function analyze_VCsaw(filename,savefolder)
%% Made by Ate Bijlsma, s4212215, a.bijlsma@neurophysiology.nl | ate.bijlsma@student.ru.nl

% analyze_VCsaw is part of the mother function 'analyze_patch' and is a function 
% used for a voltage clamp sawtooth protocol analysis

% This function retrieves the data exported from HEKA patch clamp setups
% EPC 9 and EPC 10 double. It puts the data into matlab and performs data
% analysis. That data that is used has to be collected during voltage clamp
% experiments.

% Command entry example: analyze_VCsaw('KK_MagnetoSlices_12012018_E1_5_Sawtooth_WithoutMagnet')

close all; 

load ([filename '.mat']);

files = whos ('Trace*');                            % lists the current variables, in long format form
output.file.name = filename;                         
disp (['analyze_VCsaw: ' output.file.name])            % displays what you are doing, the ST analysis

%% Organization of HEKA patch clamp data:

% Trace_X_Y_Z_W
% X: experiment number
% Y: protocol number
% Z: sweep
% W: channels.  

% Example: 
% Experiment AB8 with the protocols CC_C,CC_C,CC_C,ST100_C,ST50_C,ST10_C,VC_C in that order
% analyze_VCsaw('170926_AB_ST100_C')
% Resulting in Trace_X_Y_Z_W
% X = 8     (AB8)
% Y = 4     (ST100_C is in position 4 in the chain of protocols shown in line 31
% Z = 1-2/3 (ST100 protocol consist of 2/3 repetitions, Z will change from 1 to 3)
% W = 1-4   (Setup 1 has 4 channels (Imon2-Vmon2-Imon1-Vmon1), setup 2/3 has
%           2 channels (Vmon1-Imon1). As this experiment is performed on setup 1, we
%           will have 4 different channels. W will change from 1 to 4.

% Resulting in Trace_8_4_1_1 for the first sweep in channel 1

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

%% Channel configuration

% Depending on the setup, the traces will consist of multiple channels.
% Here we let the code know where to look for the current and voltge data.sweep{ctr}.
if round(max(output.data.values{2}(:,2)),2) == 0.0700 || round(max(output.data.values{2}(:,2)),2) == 0.0800 || round(max(output.data.values{2}(:,2)),2) == 0.1000
    Ichannel = 1;
    Vchannel = 2;
else
    Ichannel = 3;
    Vchannel = 4;
end

%% Stimulation configuration

if files(1).size(1) == 22000 || files(1).size(1) == 55000 || files(1).size(1) == 54750 % Saw 100
    stimbegin = 0.045:0.2:0.845;   
    mindistance = 0.19;    
    protocol = '100';
elseif files(1).size(1) == 42000 % Saw 200
    stimbegin = 0.045:0.4:1.645;    
    mindistance = 0.39;  
    protocol = '200';
elseif files(1).size(1) == 30000 || files(1).size(1) == 11000 || files(1).size(1) == 29750 % Saw 50
    stimbegin = 0.045:0.1:0.445;    
    mindistance = 0.09;  
    protocol = '50';
elseif files(1).size(1) == 10000 || files(1).size(1) == 9750 || files(1).size(1) == 4000 % Saw 10
    stimbegin = 0.045:0.02:0.125;    
    mindistance = 0.019; 
    protocol = '10';
end

% Sweepstartindex is calculated with the help of stimbegin which changes
% for each protocol. 
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

%% Vmon generator
% Some experiments had a missing Vmon channel. This channel can be made
% with the analyze_VC_generator function.

if isfield(output,'Vmon') == 0
    
    disp('VCgenerator inputs needed!!!')
    [output.Vmon.raw(1,:)] = analyze_VCgenerator(input('Hold potential [mV] = '),20000,1.1,0.05,0.1);

end

output.Vmon.raw(2,:) = output.Vmon.raw(1,:);
output.Vmon.raw(3,:) = output.Vmon.raw(1,:);
% This will lead to an row vector with a self made sawtooth Vmon protcol.
% The generation has to be done manually for all the experiments between
% 19122017_E1 and 21122017_E2.

%% Prelocations

for ctr = 1:max(output.exp.sweep)
    % Output.data.sweep{ctr}.plots is used for the main spikes 
    % of each sawtooth (the first sodium channel spike). These are used to set
    % some parameters and calculate
    output.data.sweep{ctr}.plots.pksindex = zeros(1,length(stimbegin));
    output.data.sweep{ctr}.plots.pksinitiation = zeros(1,length(stimbegin));
    output.data.sweep{ctr}.plots.locAfterStim = nan(1,length(stimbegin));
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

%% Collecting/analysis of the data

for ctr = 1:max(output.exp.sweep)   % loops over the analysis for the amount of sweeps (2 or 3).
    
    % The function findpeaks from the signalling toolbox is used to find
    % the main spikes of each sawotooths. This is done by implementing a
    % minimal distance between peaks. The mindistance is a bit smaller than
    % the step size of the stimulus.
    [output.data.sweep{ctr}.plots.pks(1,:),output.data.sweep{ctr}.plots.loc(1,:)] = findpeaks(-output.Imon.raw(ctr,:),output.data.sweep{ctr}.time{1},'MinPeakHeight',abs(0 * (-output.Imon.raw(ctr,1))),'MinPeakDistance',mindistance, 'MinPeakProminence', 0.1);
    [output.data.sweep{ctr}.plots.pkssaw(1,:),output.data.sweep{ctr}.plots.locsaw(1,:)] = findpeaks(output.Imon.raw(ctr,:),output.data.sweep{ctr}.time{1},'MinPeakDistance',mindistance, 'MinPeakProminence', 0.1);

    % Colleting of important parameters like the initiationtime of the
    % peak, index of pks and initiation, location of the spike relative to
    % the stimulus.
    if size(output.data.sweep{ctr}.plots.pks,2) == 5
        for i = 1:length(output.data.sweep{ctr}.plots.pks)
            output.data.sweep{ctr}.plots.pksindex(1,i)= find(-output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex(i):end) == output.data.sweep{ctr}.plots.pks(1,i),1) + output.data.sweep{ctr}.plots.sweepstartindex(i) - 1;
            output.data.sweep{ctr}.plots.pkssawindex(1,i)= find(output.Imon.raw(ctr,(output.data.sweep{ctr}.plots.pksindex(1,i):end)) == output.data.sweep{ctr}.plots.pkssaw(1,i),1) + output.data.sweep{ctr}.plots.pksindex(1,i)- 1;
            [output.data.sweep{ctr}.plots.integral(1,i),output.data.sweep{ctr}.plots.pksinitiation(1,i),output.data.sweep{ctr}.plots.pksinitiationindex(1,i),output.data.sweep{ctr}.plots.onsettime(1,i),output.data.sweep{ctr}.plots.offsettime(1,i),output.data.sweep{ctr}.plots.onsetcurrent(1,i),output.data.sweep{ctr}.plots.offsetcurrent(1,i),output.data.sweep{ctr}.plots.onsetslopecurrent(1,i),output.data.sweep{ctr}.plots.offsetslopecurrent(1,i)] = analyze_VCintegral(output.Imon.raw(ctr,output.data.sweep{ctr}.plots.sweepstartindex(i):end),output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.sweepstartindex(i):end,1)',output.data.sweep{ctr}.plots.sweepstartindex(i),protocol);
            output.data.sweep{ctr}.plots.locAfterStim(1,i)= output.data.sweep{ctr}.plots.loc(1,i)- (mean(diff(stimbegin)) * (i-1) + stimbegin(1));
        end
        
        for i = 1:length(output.data.sweep{ctr}.plots.pks)
            % Calculation of the halfwidth of the main spike.
            hh = abs((-output.Imon.raw(ctr,:)) - (output.data.sweep{ctr}.plots.pks(1,i)- ((output.data.sweep{ctr}.plots.pks(1,i)+ output.data.sweep{ctr}.plots.pksinitiation(1,i)) / 2)));
            output.data.sweep{ctr}.plots.halfwidth1(1,i)= find(hh(output.data.sweep{ctr}.plots.pksinitiationindex(1,i):output.data.sweep{ctr}.plots.pksindex(1,i)) == min(hh(output.data.sweep{ctr}.plots.pksinitiationindex(1,i):output.data.sweep{ctr}.plots.pksindex(1,i)),[],2)) + output.data.sweep{ctr}.plots.pksinitiationindex(1,i);
            output.data.sweep{ctr}.plots.halfwidth2(1,i)= (output.data.sweep{ctr}.plots.halfwidth1(1,i)- 1) + (find(output.Imon.raw(ctr,output.data.sweep{ctr}.plots.halfwidth1(1,i):end) > output.Imon.raw(ctr,output.data.sweep{ctr}.plots.halfwidth1(1,i)),1));
        end
        
        % Now all the data around every spike in the measurment is collected.
        % An if/else construction is used to place data from sweep 1 in
        % output.data.sweep{ctr}.findpeaks.sweep etc..
        for peak = 1:length(output.data.sweep{ctr}.plots.pks)
            [intE,intEl,intEw,intEp] = findpeaks(-output.Imon.raw(ctr,output.data.sweep{ctr}.plots.pksinitiationindex(1,peak):output.data.sweep{ctr}.plots.pkssawindex(1,peak)),'MinPeakProminence',0.1);
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
        
        % Loop for the timing of second peaks.
        for j = 1:length(output.data.sweep{ctr}.plots.pks)
            try
                output.data.sweep{ctr}.plots.locAfterStim2(1,j) = output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.findpeaks.sweeploc(j,2)) - (mean(diff(stimbegin)) * (j-1) + stimbegin(1));
            catch
                output.data.sweep{ctr}.plots.locAfterStim2(1,j) = NaN;
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
        
        output.data.sweep{ctr}.plots.initiation(1,:) = output.Vmon.raw(ctr,output.data.sweep{ctr}.plots.pksinitiationindex(1,:));
        output.data.sweep{ctr}.plots.halfwidth(1,:) = output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.halfwidth2(1,:)) - output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.halfwidth1(1,:));
        output.data.sweep{ctr}.plots.currentVStime(1,:) = (output.data.sweep{ctr}.plots.pks(1,:) + output.data.sweep{ctr}.plots.pksinitiation(1,:))./((output.data.sweep{ctr}.plots.loc(1,:)-(output.data.sweep{ctr}.time{1}(output.data.sweep{ctr}.plots.pksinitiationindex(1,:)))')*1000);
    else
        output.data.sweep{ctr}.plots.pks(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.loc(1,1:5) = NaN;
        output.data.sweep{ctr}.findpeaks.sweeptimediff(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.pksinitiationindex(1,1:5) = NaN;
        output.data.sweep{ctr}.findpeaks.sweepevents(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.currentVStime(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.halfwidth(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.initiation(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.locAfterStim(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.locAfterStim2(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.integral(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.onsettime(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.offsettime(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.onsetslopecurrent(1,1:5) = NaN;
        output.data.sweep{ctr}.plots.offsetslopecurrent(1,1:5) = NaN;      
     end
end

if exist([cd filesep savefolder],'dir') == 0
    mkdir(savefolder); save([savefolder '/analyzed_VCsaw_' filename], 'output')        
else
    save([savefolder '/analyzed_VCsaw_' filename], 'output')
end